#!/bin/bash

ls m* | while read id
do
  cut -f 1,4,6,10 ${id} | sed '/^@/d' >> allseq_info.sam
done