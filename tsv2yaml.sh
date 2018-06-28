#!/bin/bash
echo "Sample:"
sed '/^#/d;/^$/d' $@ |
(
while read samplename fq1 fq2 coment;do
	cat <<EOF
    ${samplename}: $comment
        fq1: $fq1
        fq2: $fq2
EOF
done
)
