docker build -t mandy/varscan varscan_dockerfile

docker run -v "/Users/myan/Documents/docker_test/bams":"/input"  mandy/varscan /usr/bin/samtools-1.9/samtools mpileup -q1 -f /input/hs37d5.fa -o /input/pileup /input/12.bam

docker exec 180c573f8890 /usr/bin/samtools-1.9/samtools mpileup -q1 -f /input/hs37d5.fa -o /input/tumor.pileup /input/Whole_T1_K1ID2_B81875.bam

docker run --rm  -v  "/Users/myan/Documents/docker_test/bams":"/input"  mandy/varscan java -jar varscan.jar somatic /input/normal.pileup /input/tumor.pileup /input/4438.varscan
