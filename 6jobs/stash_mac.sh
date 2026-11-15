#!/bin/bash
# After a sim run, moves its macro from build/ into its run folder (build/root/<the macro's /output/folder>/),
# so build/ keeps only the templates. Templates (run_0.mac, run_nophysics.mac) and failed runs are left in place.
# Usage, from build/: bash stash_mac.sh <macro> <sim exit code>
MAC=$1; STATUS=$2
[ "$STATUS" = 0 ] || exit 0
case "$(echo "$MAC" | tr 'A-Z' 'a-z')" in run_0.mac|run_nophysics.mac) exit 0 ;; esac
folder=$(awk '/^\/output\/folder /{print $2}' "$MAC" | tail -1)
[ -n "$folder" ] && [ -d "root/$folder" ] && mv "$MAC" "root/$folder/" && echo "moved $MAC to root/$folder/"
exit 0
