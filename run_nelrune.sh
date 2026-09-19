SAMPLE="ZD-4631-BcellsLaneF"
rm -rf /home/med-sal/sens05_home/NAS/NELRUNE/REAL_TEST_M39/$SAMPLE/vdj/vdj_out

target/release/nelrune-vdj \
    --bam /home/med-sal/sens05_home/NAS/NELRUNE/REAL_TEST_M39/$SAMPLE/nelrune/nelrune_out/nelrune.mapper.bam \
    --index /home/med-sal/sens05_home/NAS/NELRUNE/indexes/GRCm39_M39/mouse_GRCm39_M39.vdjidx \
    --out /home/med-sal/sens05_home/NAS/NELRUNE/REAL_TEST_M39/$SAMPLE/vdj/vdj_out \
    --exonic /home/med-sal/sens05_home/NAS/NELRUNE/REAL_TEST_M39/$SAMPLE/nelrune/nelrune_out/exonic/ \
    --threads 8 \
    --write-sequences
