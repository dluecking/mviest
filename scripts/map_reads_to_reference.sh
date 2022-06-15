# dlu 08.12.21
# needs the mapping env activate
# output: a) {READS_NAME}_vs_{REFERENCE_NAME}.mapped_READS.fastq
#         b) {READS_NAME}_vs_{REFERENCE_NAME}.coverage.tsv
# runtime couple of minutes 



READS=`realpath $1`
REFERENCE=`realpath $2`
OUTPUT_NAME=$3
WD=$4
CPUs=$5

READS_NAME=`basename $READS`
REFERENCE_NAME=`basename $REFERENCE`

echo "##########################################################"
echo "Input"
echo "reads: $READS_NAME"
echo "reference: $REFERENCE_NAME"
echo "output name: $OUTPUT_NAME"
echo "wd: $WD"

echo "------ mapping -------"

if [[ ! -e $READS ]]; then
    echo "file $READS does not exist, aborting..."
    exit
elif [[ ! -e $REFERENCE ]]; then
    echo "file $REFERENCE does not exist, aborting..."
    exit
fi

mkdir -p $WD
cd $WD

echo "1. mapping reads to reference"
bwa mem -t $CPUs -d 0 -v 1 $REFERENCE $READS > ${OUTPUT_NAME}.sam
echo "2. from sam to bam"
samtools view -@ $CPUs -bS -T $REFERENCE ${OUTPUT_NAME}.sam > ${OUTPUT_NAME}.bam
echo "3. sorting"
samtools sort -@ $CPUs -T tmp ${OUTPUT_NAME}.bam > ${OUTPUT_NAME}.sorted.bam
echo "4. calculating depth"
samtools depth ${OUTPUT_NAME}.sorted.bam > ${OUTPUT_NAME}.sorted.bam.coverage
echo "5. extracting mapped reads"
bam2fastq --no-unaligned --force --strict -o ${OUTPUT_NAME}.mapped_reads.fastq ${OUTPUT_NAME}.sorted.bam
echo "6. zipping reads and cleaning up"
pigz ${OUTPUT_NAME}.mapped_reads.fastq
rm ${OUTPUT_NAME}.sam
rm ${OUTPUT_NAME}.bam
echo "mapping done..."
