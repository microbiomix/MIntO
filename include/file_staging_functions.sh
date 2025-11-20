# ---------------------------------------------------------------------------
# file_staging_functions.sh
#
# Bash function library for safe, concurrent staging of large reference files
# from a shared filesystem (e.g. NFS) onto node-local storage in a cluster.
#
# This file is *not* meant to be executed directly. Instead, source it:
#
#     source "{minto_dir}/include/file_staging_functions.sh"
#
# It provides the following functions:
#
#   stage_file_in LOCAL_DIR REMOTE_FILE
#       Safely stage a single file from a shared path to a node-local path.
#       Uses per-file flock-based locking to ensure that:
#           - only one rsync occurs per file per node,
#           - subsequent jobs reuse the staged local copy,
#           - concurrent jobs do not overload the NFS,
#           - lock duration is short and limited to the rsync window.
#
#   stage_multiple_files_in LOCAL_DIR REMOTE_FILE1 [REMOTE_FILE2 ...]
#       Convenience wrapper around stage_file_in for staging several files.
#       Includes jitter to reduce simultaneous lock contention.
#
# Locks:
#   Locks are implemented via flock() on file descriptors, not lockfile(-like)
#   staleness conventions. This ensures:
#       - automatic lock release on process exit,
#       - no stale lockfiles,
#       - race-free mutual exclusion even under heavy cluster load.
#
# Assumptions:
#   - LOCAL_DIR is node-local storage (e.g. /scratch, /tmp, /local).
#   - REMOTE_FILE is accessible via shared FS (NFS / Lustre / GPFS).
#   - rsync is available on execution nodes.
#
# ---------------------------------------------------------------------------

#################################
# Disable direct execution
#################################

if [[ "${BASH_SOURCE[0]}" == "${0}" ]]; then
    echo "ERROR: This file must be sourced, not executed." >&2
    exit 1
fi

#################################
# Function definitions
#################################

# Stage one file to the local location using rsync, with flock-based locking.
# Usage: stage_file_in LOCAL_DIR LOCK_ID REMOTE_FILE
stage_file_in() {
    local local_dir=$1
    local lock_id=$2
    local remote_file=$3

    if [[ -z "$local_dir" || -z "$lock_id" || -z "$remote_file" ]]; then
        echo "Usage: stage_file_in LOCAL_DIR LOCK_ID REMOTE_FILE" >&2
        return 2
    fi

    # Make sure the target directory exists
    mkdir -p -- "$local_dir" || {
        echo "Cannot create directory: $local_dir" >&2
        return 1
    }

    local file local_copy remote_lock local_lock
    file=$(basename -- "$remote_file")
    local_copy="${local_dir%/}/$file"

    # Lock filenames (can live on NFS or local)
    remote_lock="${remote_file}.lock.${lock_id}"
    local_lock="${local_copy}.lock"

    # Acquire an exclusive lock on the remote lockfile
    exec 9> "$remote_lock" || {
        echo "Cannot open lock file $remote_lock" >&2
        return 1
    }
    # Wait up to 3600s to acquire the lock
    if ! flock -w 3600 9; then
        echo "Could not acquire lock on $remote_lock within 3600 seconds" >&2
        exec 9>&-
        return 1
    fi

    # Acquire an exclusive lock on the local copy lockfile
    exec 8> "$local_lock" || {
        echo "Cannot open lock file $local_lock" >&2
        exec 9>&-
        return 1
    }
    if ! flock -w 3600 8; then
        echo "Could not acquire lock on $local_lock within 3600 seconds" >&2
        exec 8>&-
        exec 9>&-
        return 1
    fi

    # At this point we hold both locks (FDs 9 and 8).
    # But ONLY for the duration of the rsync loop.

    local attempts=0
    local rsync_rc=0

    # We set max_attempts to 1, as it is designed for use from within Snakemake
    # where "fail early" is the best approach. There might already be a
    # "--restart-times=N" option given to Snakemake, which might make this
    # even slower to catch an error when max_attempts > 1.
    while (( attempts < 1 )); do
        if rsync --itemize-changes -a -- "$remote_file" "$local_copy"; then
            rsync_rc=0
            break
        fi
        rsync_rc=$?
        attempts=$(( attempts + 1 ))
        echo "rsync failed for $remote_file -> $local_copy (attempt $attempts, rc=$rsync_rc), retrying..." >&2
        sleep 1
    done

    # release the locks before returning
    exec 8>&-
    exec 9>&-

    if (( attempts == 10 )); then
        echo "Staging $remote_file to $local_copy failed after $attempts attempts (last rc=$rsync_rc)." >&2
        echo "Please copy it manually or investigate network/storage issues." >&2
        return 1
    fi

    # No need to rm lock files. The locks are released automatically when
    # the shell closes FDs 8 and 9 (i.e. when this function/parent shell exits).

    return 0
}

# Stage a list of files to the local location
# Usage: stage_multiple_files_in LOCAL_DIR REMOTE_FILE1 [REMOTE_FILE2 ...]
stage_multiple_files_in() {
    local local_location=$1

    # Require at least 2 arguments: local dir, and at least one remote file
    if [[ -z "$local_location" || $# -lt 2 ]]; then
        echo "Usage: stage_multiple_files_in LOCAL_DIR LOCK_ID REMOTE_FILE1 [REMOTE_FILE2 ...]" >&2
        return 2
    fi

    # Skip the first positional argument (local_location)
    shift 1

    # Make local dir if it doesn't exist
    mkdir -p -- "$local_location" || {
        echo "Cannot create directory: $local_location" >&2
        return 1
    }

    # Jitter: wait for a random amount of seconds (1–60) to avoid stampedes
    #sleep $((RANDOM % 60 + 1))

    # Get one of five host-specific locks randomly
    # Determine stable lock id based on host name
    local remote
    local host=${HOSTNAME:-$(hostname)}
    local md5=$(printf '%s' "$host" | md5sum | awk '{print $1}')
    local lock_id=$(( 0x${md5:0:4} % 5 + 1 ))

    # Stage files one by one
    for remote in "$@"; do
        echo "Staging: $remote -> $local_location" >&2
        # Let stage_file_in handle its own mkdir/locking/rsync logic
        if ! stage_file_in "$local_location" "$lock_id" "$remote"; then
            echo "Failed to stage $remote to $local_location" >&2
            return 1
        fi
    done

    return 0
}

stage_file_out() {
	local local_file=$1;
	local remote_dir=$2;
	if [ -z "$remote_dir" ]; then
		echo "Usage: stage_file_out local-file remote-dir" 1>&2;
		exit 2;
	fi;

	local file=$(basename $local_file);
	local remote_copy="$remote_dir/$file";

	local attempts=0;
	local copy=1;
	while [ $attempts -lt 10 -a $copy -eq 1 ]; do
		rsync -q -ptg $local_file $remote_copy;
		if [ $? -ne 0 ] ; then  # rsync failed!
			copy=1;
			attempts=$(expr $attempts + 1);
		else
			copy=0;
		fi;
	done
	if [ $attempts -eq 10 ]; then
		echo "Staging out $local_file to $remote_dir/ failed after $attempts attempts! Please copy it yourself by running:" 1>&2;
		echo "rsync -a $(hostname):$local_file $remote_dir/" 1>&2;
	else
		rm -f $local_file;
	fi
}
