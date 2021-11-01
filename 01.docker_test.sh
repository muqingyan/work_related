docker build -t mandy/varscan varscan_dockerfile

docker run -v "/Users/myan/Documents/Natera/docker_test/bams":"/input"  mandy/varscan /usr/bin/samtools-1.9/samtools mpileup -q1 -f /input/hs37d5.fa -o /input/C032_4438.normal.pileup /input/12050926.3-2-BXE_C032_4438_034423_PB_Whole_C1_K1ID2_B81912.bam

docker exec 180c573f8890 /usr/bin/samtools-1.9/samtools mpileup -q1 -f /input/hs37d5.fa -o /input/C032_4438.tumor.pileup /input/12513698.1-2-FXB_C032_4438_035042_LV_Whole_T1_K1ID2_B81875.bam

docker run --rm  -v  "/Users/myan/Documents/Natera/docker_test/bams":"/input"  mandy/varscan java -jar varscan.jar somatic /input/C032_4438.normal.pileup /input/C032_4438.tumor.pileup /input/C032_4438.varscan
