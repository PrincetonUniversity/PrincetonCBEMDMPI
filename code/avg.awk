#!/usr/bin/env awk -f
# awk -v field=4 -f avg.awk cool > file   will tell you the average value

# need to run ./andersen ... > file 
# then need to grep lines with "Instant" to cool
# then run this 

{ tot += $field; count++ }
END { print tot/count }
