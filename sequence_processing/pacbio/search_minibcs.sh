#!/bin/bash

# Script to search a file for a set of substrings and output the count of each

usage() {
    echo "Usage: $0 <file> <substring1> [substring2] [substring3] ..."
    echo ""
    echo "Searches <file> for each substring and outputs the count of occurrences."
    echo ""
    echo "Example:"
    echo "  $0 myfile.txt 'hello' 'world' 'foo'"
    exit 1
}

# Check for minimum arguments
if [[ $# -lt 2 ]]; then
    usage
fi

file="$1"
shift

# Check if file exists and is readable
if [[ ! -f "$file" ]]; then
    echo "Error: File '$file' not found." >&2
    exit 1
fi

if [[ ! -r "$file" ]]; then
    echo "Error: File '$file' is not readable." >&2
    exit 1
fi

# Search for each substring and count occurrences
for substring in "$@"; do
    # -o outputs each match on its own line
    # -F treats pattern as fixed string (not regex)
    count=$(grep -o -F -- "$substring" "$file" 2>/dev/null | wc -l)
    echo "$substring: $count"
done
