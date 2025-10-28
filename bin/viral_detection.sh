#!/bin/bash

# Check if the input file is provided
if [[ $# -ne 1 ]]; then
    echo "Usage: $0 <input_file>"
    exit 1
fi

# Read the input file
INPUT_FILE="$1"

if [[ ! -f "$INPUT_FILE" ]]; then
    echo "Error: File '$INPUT_FILE' not found."
    exit 1
fi

# Parse the input file
while IFS="=" read -r key value; do
    case "$key" in
	RNA)
            RNA="$value"
            ;;
	DNA)
            DNA="$value"
            ;;
        OUTPUT_FOLDER)
            OUTPUT_FOLDER="$value"
            ;;
        GEX_FASTQ_FOLDER)
            GEX_FASTQ_FOLDER="$value"
            ;;
	ATAC_FASTQ_FOLDER)
            ATAC_FASTQ_FOLDER="$value"
            ;;
        GEX_WHITE)
            GEX_WHITE="$value"
            ;;
        ATAC_WHITE)
            ATAC_WHITE="$value"
            ;;
        CORES)
            CORES="$value"
            ;;
        VIRAL_INDEX_FOLDER_GEX)
            VIRAL_INDEX_FOLDER_GEX="$value"
            ;;
        VIRAL_INDEX_FOLDER_ATAC)
            VIRAL_INDEX_FOLDER_ATAC="$value"
            ;;
	VIRAL_INDEX_NAME)
            VIRAL_INDEX_NAME="$value"
            ;;
        VIRAL_FASTA)
            VIRAL_FASTA="$value"
            ;;
        ADJUSTMENT_FACTOR)
            ADJUSTMENT_FACTOR="$value"
            ;;
	GEX_READ_GREP)
	    GEX_READ_GREP="$value"
	    ;;
	GEX_CB_GREP)
	    GEX_CB_GREP="$value"
	    ;;
        ATAC_READ1_GREP)
            ATAC_READ1_GREP="$value"
            ;;
        ATAC_READ2_GREP)
            ATAC_READ2_GREP="$value"
            ;;
        ATAC_CB_GREP)
            ATAC_CB_GREP="$value"
            ;;
	ATAC_SIGNIFIER)
            ATAC_SIGNIFIER="$value"
            ;;
        *)
            echo "Warning: Unknown parameter '$key' in the input file."
            ;;
    esac
done < "$INPUT_FILE"

# Create an array to store file names
declare -a GEX_FILE_NAMES

# Create an array to store file names
declare -a ATAC_FILE_NAMES

running_jobs=0

if [[ ! -d "$OUTPUT_FOLDER" ]]; then
    mkdir $OUTPUT_FOLDER
fi

# Check for ambiguous bases in the viral FASTA
if grep -q '[^ACGTNacgtn>[:space:]]' "$VIRAL_FASTA"; then
    echo "⚠️  Ambiguous bases detected in $VIRAL_FASTA — creating cleaned FASTA..."

    CLEAN_FASTA="${VIRAL_FASTA%.fa*}_clean.fa"

    # Replace any non-ACGTN bases with 'N'
    awk '/^>/ {print; next} {gsub(/[^ACGTNacgtn]/,"N"); print}' "$VIRAL_FASTA" > "$CLEAN_FASTA"

    # Index the cleaned FASTA
    samtools faidx "$CLEAN_FASTA"

    # Update the variable so downstream commands use the cleaned reference
    VIRAL_FASTA="$CLEAN_FASTA"

    echo "✅ Cleaned FASTA created and indexed: $CLEAN_FASTA"
else
    echo "✅ No ambiguous bases detected in $VIRAL_FASTA — using original reference."
fi

if [ "$RNA" = "TRUE" ]; then

    # Check if viral index folder exists
    if [[ ! -d "$VIRAL_INDEX_FOLDER_GEX" ]]; then
        STAR --runThreadN $CORES \
	     --runMode genomeGenerate \
	     --genomeDir $VIRAL_INDEX_FOLDER_GEX \
	     --genomeFastaFiles $VIRAL_FASTA \
	     --genomeSAindexNbases $ADJUSTMENT_FACTOR
    fi

    IFS=',' read -r -a GEX_FOLDERS <<< "$GEX_FASTQ_FOLDER"
    GEX_FILE_NAMES=()

    # Find files matching the naming convention and store them in array
    for folder in "${GEX_FOLDERS[@]}"; do
        while IFS= read -r -d $'\0' file; do
            FILE_NAME=$(basename "$file")
            NAME_PART=$(echo "$FILE_NAME" | awk -F '_L[0-9]*' '{print $1}')

            # Skip if the name part contains "ATAC"
            if [[ ! $NAME_PART =~ $ATAC_SIGNIFIER ]]; then
                GEX_FILE_NAMES+=("$NAME_PART|$FILE_NAME|$folder")
            fi
        done < <(find -L "$folder" -type f -name "$GEX_READ_GREP" -print0)
    done

    # Calculate the largest possible number of parallel runs given a minimum of 8 cores per run
    max_samples=$((CORES / 8))
    max_samples=$((max_samples < ${#GEX_FILE_NAMES[@]} ? max_samples : ${#GEX_FILE_NAMES[@]}))

    # Determine number of cores and memory per run
    cores_per_sample=$((CORES / max_samples))

    # Run STAR for each file
    for file_info in "${GEX_FILE_NAMES[@]}"; do
        IFS="|" read -r NAME_PART FILE_NAME FOLDER <<< "$file_info"
        if [[ ! -d "$OUTPUT_FOLDER/$NAME_PART" ]]; then
            mkdir $OUTPUT_FOLDER/$NAME_PART
        fi
        if [ ! -e "$OUTPUT_FOLDER/$NAME_PART/${NAME_PART}Log.out" ]; then
            STAR --genomeDir $VIRAL_INDEX_FOLDER_GEX \
		 --runThreadN $cores_per_sample \
		 --readFilesIn $FOLDER/$FILE_NAME \
		 --readFilesCommand zcat \
		 --outFilterMultimapNmax 100000 \
		 --outFilterScoreMin 30 \
		 --outSAMmultNmax -1 \
		 --outMultimapperOrder Random \
		 --outFileNamePrefix $OUTPUT_FOLDER/$NAME_PART/$NAME_PART \
		 --outSAMtype BAM SortedByCoordinate &
        else
	    continue
        fi
        ((running_jobs++))

        # Check to see if as many jobs as possible are running
        if [ "$running_jobs" -eq "$max_samples" ]; then
            wait -n
            ((running_jobs--))
        fi
    done

    wait

    # Run bcftools for each file to create a consensus fasta
    for file_info in "${GEX_FILE_NAMES[@]}"; do
        IFS="|" read -r NAME_PART FILE_NAME FOLDER <<< "$file_info"
        if [[ ! -e "$OUTPUT_FOLDER/$NAME_PART/${NAME_PART}_variants.norm.vcf.gz" ]]; then
            samtools index $OUTPUT_FOLDER/$NAME_PART/${NAME_PART}Aligned.sortedByCoord.out.bam
            bcftools mpileup -f $VIRAL_FASTA $OUTPUT_FOLDER/$NAME_PART/${NAME_PART}Aligned.sortedByCoord.out.bam -Ou |
            bcftools call -mv -Ou |
            bcftools norm -f $VIRAL_FASTA -Oz -o $OUTPUT_FOLDER/$NAME_PART/${NAME_PART}_variants.norm.vcf.gz
	    bcftools index $OUTPUT_FOLDER/$NAME_PART/${NAME_PART}_variants.norm.vcf.gz
        fi
    done

fi

if [ "$DNA" = "TRUE" ]; then

    # Check if viral index folder exists
    if [[ ! -d "$VIRAL_INDEX_FOLDER_ATAC" ]]; then
        mkdir $VIRAL_INDEX_FOLDER_ATAC
        bowtie2-build $VIRAL_FASTA "$VIRAL_INDEX_FOLDER_ATAC/$VIRAL_INDEX_NAME"
    fi

    IFS=',' read -r -a ATAC_FOLDERS <<< "$ATAC_FASTQ_FOLDER"
    ATAC_FILE_NAMES=()

    # Find files matching the naming convention and store them in array
    for folder in "${ATAC_FOLDERS[@]}"; do
        while IFS= read -r -d $'\0' file; do
            FILE_NAME=$(basename "$file")

            # Extract the name part
            NAME_PART=$(echo "$FILE_NAME" | awk -F '_L[0-9]*' '{print $1}')

            # Check if name part contains "ATAC"
            if [[ $NAME_PART =~ $ATAC_SIGNIFIER ]]; then
                ATAC_FILE_NAMES+=("$NAME_PART|$FILE_NAME|$folder")
            fi
        done < <(find -L "$folder" -type f -name "$ATAC_READ1_GREP" -print0)
    done

    # Calculate the largest possible number of parallel runs given a minimum of 8 cores per run
    max_samples=$((CORES / 8))
    max_samples=$((max_samples < ${#ATAC_FILE_NAMES[@]} ? max_samples : ${#ATAC_FILE_NAMES[@]}))

    # Determine number of cores and memory per run
    cores_per_sample=$((CORES / max_samples))

    running_jobs=0

    # Run bowtie2 for each file
    for file_info in "${ATAC_FILE_NAMES[@]}"; do
        IFS="|" read -r NAME_PART FILE_NAME FOLDER <<< "$file_info"
        if [[ ! -d "$OUTPUT_FOLDER/$NAME_PART" ]]; then
            mkdir $OUTPUT_FOLDER/$NAME_PART
        fi
        if [ ! -e "$OUTPUT_FOLDER/$NAME_PART/output.sorted.bam" ]; then
            # Find related FASTQ file
	    for folder in "${ATAC_FOLDERS[@]}"; do
    		RELATED_FILE=$(find -L "$folder" -type f -name "${NAME_PART}${ATAC_READ2_GREP}" -print -quit)
    		if [[ -n "$RELATED_FILE" ]]; then
        	    break
	        fi
	    done
            bowtie2 --no-discordant \
		    --local -a\
		    -x "$VIRAL_INDEX_FOLDER_ATAC/$VIRAL_INDEX_NAME" \
		    -1 $FOLDER/$FILE_NAME \
		    -2 $RELATED_FILE \
		    -N 1 --very-sensitive \
		    -p $cores_per_sample --no-unal |
		    samtools view -bS |
		    samtools sort -o "$OUTPUT_FOLDER/$NAME_PART/output.sorted.bam" &
        else
            continue
        fi
        ((running_jobs++))

        # Check to see if as many jobs as possible are running
        if [ "$running_jobs" -eq "$max_samples" ]; then
            wait -n
            ((running_jobs--))
        fi
    done

    wait

    # Run bcftools for each file with to create a consensus fasta
    for file_info in "${ATAC_FILE_NAMES[@]}"; do
        IFS="|" read -r NAME_PART FILE_NAME FOLDER <<< "$file_info"
        if [[ ! -e "$OUTPUT_FOLDER/$NAME_PART/${NAME_PART}_variants.norm.vcf.gz" ]]; then
            samtools index $OUTPUT_FOLDER/$NAME_PART/output.sorted.bam
            bcftools mpileup -f $VIRAL_FASTA $OUTPUT_FOLDER/$NAME_PART/output.sorted.bam -Ou |
            bcftools call -mv -Ou |
            bcftools norm -f $VIRAL_FASTA -Oz -o $OUTPUT_FOLDER/$NAME_PART/${NAME_PART}_variants.norm.vcf.gz
	    bcftools index $OUTPUT_FOLDER/$NAME_PART/${NAME_PART}_variants.norm.vcf.gz
        fi
    done

fi

# Collect prefixes before first underscore
prefixes=()
for f in "${GEX_FILE_NAMES[@]}" "${ATAC_FILE_NAMES[@]}"; do
    prefixes+=( "${f%%|*}" )
    prefixes[-1]="${prefixes[-1]%%_*}"
done

# Get unique prefixes
readarray -t unique_prefixes < <(printf "%s\n" "${prefixes[@]}" | sort -u)

# Calculate the largest possible number of parallel runs given a minimum of 8 cores per run
max_samples=$((CORES / 8))
max_samples=$((max_samples < ${#unique_prefixes[@]} ? max_samples : ${#unique_prefixes[@]}))

# Determine number of cores and memory per run
cores_per_sample=$((CORES / max_samples))

# Define the function that processes one prefix
process_prefix() {
    local prefix="$1"
    local files_to_merge=()
    local has_gex=false
    local has_atac=false
    local GEX_FILE_NAMES=($GEX_FILE_NAMES_STR)
    local ATAC_FILE_NAMES=($ATAC_FILE_NAMES_STR)

    # Collect matching files and track type
    for f in "${GEX_FILE_NAMES[@]}"; do
        if [[ "$f" == "${prefix}"_* ]]; then
	    IFS="|" read -r NAME_PART FILE_NAME FOLDER <<< "$f"
            files_to_merge+=("$OUTPUT_FOLDER/$NAME_PART/${NAME_PART}_variants.norm.vcf.gz")
            has_gex=true
        fi
    done
    for f in "${ATAC_FILE_NAMES[@]}"; do
        if [[ "$f" == "${prefix}"_* ]]; then
	    IFS="|" read -r NAME_PART FILE_NAME FOLDER <<< "$f"
            files_to_merge+=("$OUTPUT_FOLDER/$NAME_PART/${NAME_PART}_variants.norm.vcf.gz")
            has_atac=true
        fi
    done

    if [[ ${#files_to_merge[@]} -ne 0 ]]; then
        mkdir -p "$OUTPUT_FOLDER/$prefix"

	consensus_fasta="$OUTPUT_FOLDER/$prefix/${prefix}_viral_genomes.fasta"
	filtered_vcf="$OUTPUT_FOLDER/$prefix/${prefix}_variants.norm.vcf.gz"

	if [[ ${#files_to_merge[@]} -gt 1 ]]; then
            bcftools merge "${files_to_merge[@]}" -Ou | \
                bcftools norm -f "$VIRAL_FASTA" -m -any -Ou | \
                bcftools filter -e 'INFO/DP4[2]+INFO/DP4[3] <= INFO/DP4[0]+INFO/DP4[1]' \
                    -Oz -o "$filtered_vcf"
	else
	    bcftools norm -f "$VIRAL_FASTA" -m -any -Ou "${files_to_merge[0]}" | \
	        bcftools filter -e 'INFO/DP4[2]+INFO/DP4[3] <= INFO/DP4[0]+INFO/DP4[1]' \
        	-Oz -o "$filtered_vcf"
	fi


        bcftools index "$filtered_vcf"

        bcftools consensus -f "$VIRAL_FASTA" "$filtered_vcf" > "$consensus_fasta"

        # If this prefix was in GEX or ATAC arrays, build genome indices accordingly
        if $has_gex; then
	    if [[ ! -d "$OUTPUT_FOLDER/$prefix/STAR_index" ]]; then
                STAR --runThreadN $cores_per_sample \
		     --runMode genomeGenerate \
                     --genomeDir "$OUTPUT_FOLDER/$prefix/STAR_index" \
                     --genomeFastaFiles "$consensus_fasta" \
		     --genomeSAindexNbases $ADJUSTMENT_FACTOR
	    fi
        fi

        if $has_atac; then
            if [[ ! -d "$OUTPUT_FOLDER/$prefix/bowtie2_index" ]]; then
                mkdir "$OUTPUT_FOLDER/$prefix/bowtie2_index"
                bowtie2-build --threads $cores_per_sample \
		    	      "$consensus_fasta" "$OUTPUT_FOLDER/$prefix/bowtie2_index/${prefix}_viral"
	    fi
        fi
    fi
}

export -f process_prefix
export OUTPUT_FOLDER VIRAL_FASTA ADJUSTMENT_FACTOR cores_per_sample GEX_FILE_NAMES_STR="${GEX_FILE_NAMES[*]}" ATAC_FILE_NAMES_STR="${ATAC_FILE_NAMES[*]}"

# Use parallel to process prefixes in parallel, limiting to MAX_PARALLEL_JOBS
parallel -j $max_samples process_prefix ::: "${unique_prefixes[@]}"

if [ "$RNA" = "TRUE" ]; then

    # Calculate the largest possible number of parallel runs given a minimum of 8 cores per run
    max_samples=$((CORES / 8))
    max_samples=$((max_samples < ${#GEX_FILE_NAMES[@]} ? max_samples : ${#GEX_FILE_NAMES[@]}))

    # Determine number of cores and memory per run
    cores_per_sample=$((CORES / max_samples))

    running_jobs=0

    # Run STAR again for each file with variants
    for file_info in "${GEX_FILE_NAMES[@]}"; do
        IFS="|" read -r NAME_PART FILE_NAME FOLDER <<< "$file_info"
	prefix="$(basename "$NAME_PART" | cut -d'_' -f1)"
        if [ ! -e "$OUTPUT_FOLDER/$NAME_PART/${NAME_PART}_VCFLog.out" ]; then
	    STAR --genomeDir "$OUTPUT_FOLDER/$prefix/STAR_index" \
		 --runThreadN $cores_per_sample \
		 --readFilesIn $FOLDER/$FILE_NAME \
		 --readFilesCommand zcat \
		 --outFilterMultimapNmax 100000 \
		 --outFilterScoreMin 30 \
		 --outSAMmultNmax 1 \
		 --outFileNamePrefix $OUTPUT_FOLDER/$NAME_PART/${NAME_PART}_VCF \
		 --outSAMtype BAM Unsorted &
        else
           continue
        fi
        ((running_jobs++))

        # Check to see if as many jobs as possible are running
        if [ "$running_jobs" -eq "$max_samples" ]; then
            wait -n
            ((running_jobs--))
        fi
    done

    wait

    # Process SAM files
    process_file_gex() {
        file_info=$1
        IFS="|" read -r NAME_PART _ <<< "$file_info"
	for folder in "${GEX_FOLDERS[@]}"; do
    	    RELATED_FILE=$(find -L "$folder" -type f -name "${NAME_PART}${GEX_CB_GREP}" -print -quit)
    	    if [[ -n "$RELATED_FILE" ]]; then
        	break
    	    fi
	done
        samtools view $OUTPUT_FOLDER/$NAME_PART/${NAME_PART}_VCFAligned.out.bam |
        grep -v "GGGGGGGGGGGGGGGGGGGG" |
        awk '{print $1}' |
        sort -u |
        sed 's/\/[12]$//' |
        seqtk subseq $RELATED_FILE - |

        # Main awk script
        awk -v gex_white="$GEX_WHITE" '
        function hamming_distance(str1, str2) {
            # Calculate Hamming distance between two strings
            dist = 0
            for (i = 1; i <= length(str1); i++) {
                if (substr(str1, i, 1) != substr(str2, i, 1)) {
                    dist++
                    if (dist > 1) return dist # Early exit if distance > 1
                }
            }
            return dist
        }

        BEGIN {
            # Read GEX_WHITE into an array
            while ((getline line < gex_white) > 0) {
                gex_map[line] = 1
            }
            close(gex_white)
        }

        NR % 4 == 2 {
            # Check if the whole line has been seen before
            if (!line_seen[$0]++) {
                # Extract the key from the first 16 characters
                key = substr($0, 1, 16)
                seen[key]++
            }
        }

        END {
            # Sort keys in descending order of counts
            m = asorti(seen, sorted_keys, "@val_num_desc")

            for (j = 1; j <= m; j++) {
                key = sorted_keys[j]
                count = seen[key]

                # Check for an exact match in preloaded atac_map
                if (key in gex_map) {
                    final_map[key] += count
                    continue
                }

                approx_found = 0
                match_key = ""

                # Check previously processed keys in final_map for approximate matches
                for (prev_key in final_map) {
                    if (final_map[prev_key] > 1 && hamming_distance(key, prev_key) == 1) {
                        approx_found++
                        if (approx_found > 1) break
                        match_key = prev_key
                    }
                }

                if (approx_found > 1) continue

                if (approx_found == 1) {
                    final_map[match_key] += count
                    continue
                }

                # Check for approximate matches in atac_map
                for (gex_key in gex_map) {
                    if (hamming_distance(gex_key, key) == 1) {
                        approx_found++
                        if (approx_found > 1) break
                        match_key = gex_key
                    }
                }

                if (approx_found == 1) {
                    final_map[match_key] += count
                }
            }

            # Print results
            for (key in final_map) {
                print key ", " final_map[key]
            }
        }' | \
        sort -t, -k2,2nr > $OUTPUT_FOLDER/$NAME_PART/${NAME_PART}_Wei_unfiltered.csv
    }

    export -f process_file_gex
    export OUTPUT_FOLDER GEX_CB_GREP GEX_WHITE GEX_FOLDERS

    # Use parallel to process files in parallel, limiting to MAX_PARALLEL_JOBS
    parallel -j $max_samples process_file_gex ::: "${GEX_FILE_NAMES[@]}"
fi

if [ "$DNA" = "TRUE" ]; then
    # Calculate the largest possible number of parallel runs given a minimum of 8 cores per run
    max_samples=$((CORES / 8))
    max_samples=$((max_samples < ${#ATAC_FILE_NAMES[@]} ? max_samples : ${#ATAC_FILE_NAMES[@]}))

    # Determine number of cores and memory per run
    cores_per_sample=$((CORES / max_samples))

    running_jobs=0

    # Run bowtie2 for each file
    for file_info in "${ATAC_FILE_NAMES[@]}"; do
        IFS="|" read -r NAME_PART FILE_NAME FOLDER <<< "$file_info"
	prefix="$(basename "$NAME_PART" | cut -d'_' -f1)"
        if [ ! -e "$OUTPUT_FOLDER/$NAME_PART/output_vcf.bam" ]; then
            # Find related FASTQ file
	    for folder in "${ATAC_FOLDERS[@]}"; do
    		RELATED_FILE=$(find -L "$folder" -type f -name "${NAME_PART}${ATAC_READ2_GREP}" -print -quit)
    		if [[ -n "$RELATED_FILE" ]]; then
        	    break
    		fi
	    done
            bowtie2 --no-discordant \
                    --local -k 1 \
                    -x "$OUTPUT_FOLDER/$prefix/bowtie2_index/${prefix}_viral" \
                    -1 $FOLDER/$FILE_NAME \
                    -2 $RELATED_FILE \
                    -N 1 --very-sensitive \
                    -p $cores_per_sample --no-unal |
                    samtools view -bS -o "$OUTPUT_FOLDER/$NAME_PART/output_vcf.bam" &
        else
            continue
        fi
        ((running_jobs++))

        # Check to see if as many jobs as possible are running
        if [ "$running_jobs" -eq "$max_samples" ]; then
            wait -n
            ((running_jobs--))
        fi
    done

    wait

    # Function to process each file
    process_file_atac() {
        file_info=$1
        IFS="|" read -r NAME_PART _ <<< "$file_info"

        # Extract unique read names, sequences, count the occurrences of the first 16 characters, and format the output
        samtools view -h $OUTPUT_FOLDER/${NAME_PART}/output_vcf.bam | \
        awk '$0 !~ /^@/ && $0 !~ /GGGGGGGGGGGGGGG/ && $0 !~ /TTTTTTTTTTTTTTT/ && $0 !~ /AAAAAAAAAAAAAAA/ && $0 !~ /CCCCCCCCCCCCCCC/ {if (!seen[$1]++) print $1}' > $OUTPUT_FOLDER/${NAME_PART}/match_labels.txt
	for folder in "${ATAC_FOLDERS[@]}"; do
    	    RELATED_FILE=$(find -L "$folder" -type f -name "${NAME_PART}${ATAC_CB_GREP}" -print -quit)
    	    if [[ -n "$RELATED_FILE" ]]; then
        	break
    	    fi
	done
        seqtk subseq $RELATED_FILE $OUTPUT_FOLDER/${NAME_PART}/match_labels.txt > $OUTPUT_FOLDER/${NAME_PART}/match_R2.txt

	for folder in "${ATAC_FOLDERS[@]}"; do
    	    RELATED_FILE=$(find -L "$folder" -type f -name "${NAME_PART}${ATAC_READ2_GREP}" -print -quit)
    	    if [[ -n "$RELATED_FILE" ]]; then
        	break
    	    fi
	done
        seqtk subseq $RELATED_FILE $OUTPUT_FOLDER/${NAME_PART}/match_labels.txt > $OUTPUT_FOLDER/${NAME_PART}/match_R1.txt
        seqtk subseq $FOLDER/$FILE_NAME $OUTPUT_FOLDER/${NAME_PART}/match_labels.txt > $OUTPUT_FOLDER/${NAME_PART}/match_R3.txt

	files=("$OUTPUT_FOLDER/${NAME_PART}/match_R3.txt" "$OUTPUT_FOLDER/${NAME_PART}/match_R2.txt" "$OUTPUT_FOLDER/${NAME_PART}/match_R1.txt")

        # Output files
        output_files=("$OUTPUT_FOLDER/${NAME_PART}/filtered_match_R3_viral_unique.txt" "$OUTPUT_FOLDER/${NAME_PART}/filtered_match_R2_viral_unique.txt" "$OUTPUT_FOLDER/${NAME_PART}/filtered_match_R1_viral_unique.txt")

        # Temporary arrays to hold triplet lines
        declare -a lines1 lines2 lines3

        # Set to store unique triplets
        declare -A seen_triplets

        # Ensure output files are empty
        for output_file in "${output_files[@]}"; do
            > "$output_file"
        done

        # Open all files in parallel and process them
        exec 3<"${files[0]}" 4<"${files[1]}" 5<"${files[2]}"

        while true; do
            # Read 4-line blocks from each file
            read -r line1_1 <&3 || break
            read -r line2_1 <&3
            read -r line3_1 <&3
            read -r line4_1 <&3

            read -r line1_2 <&4 || break
            read -r line2_2 <&4
            read -r line3_2 <&4
            read -r line4_2 <&4

            read -r line1_3 <&5 || break
            read -r line2_3 <&5
            read -r line3_3 <&5
            read -r line4_3 <&5

            # Create two partial triplets
            partial_triplet1="$line2_1|${line2_2: -16}"
            partial_triplet2="${line2_2: -16}|$line2_3"

            # Check if either partial triplet has been seen before
            if [[ -z "${seen_triplets[$partial_triplet1]}" && -z "${seen_triplets[$partial_triplet2]}" && ! $partial_triplet1 =~ (G{15}|C{15}|T{15}|A{15}) && ! $partial_triplet2 =~ (G{15}|C{15}|T{15}|A{15}) ]]; then
                # Mark both partial triplets as seen
                seen_triplets["$partial_triplet1"]=1
                seen_triplets["$partial_triplet2"]=1

                # Append the blocks to the output files
                {
                    echo "$line1_1"
                    echo "$line2_1"
                    echo "$line3_1"
                    echo "$line4_1"
                } >> "${output_files[0]}"

                {
                    echo "$line1_2"
                    echo "$line2_2"
                    echo "$line3_2"
                    echo "$line4_2"
                } >> "${output_files[1]}"

                {
                    echo "$line1_3"
                    echo "$line2_3"
                    echo "$line3_3"
                    echo "$line4_3"
                } >> "${output_files[2]}"
            fi
        done

        # Close file descriptors
        exec 3<&- 4<&- 5<&-

        cat "${output_files[1]}" | \
        # Main awk script
        awk -v atac_white="$ATAC_WHITE" -v gex_white="$GEX_WHITE" '
        function hamming_distance(str1, str2) {
            # Calculate Hamming distance between two strings
            dist = 0
            for (i = 1; i <= length(str1); i++) {
                if (substr(str1, i, 1) != substr(str2, i, 1)) {
                    dist++
                    if (dist > 1) return dist # Early exit if distance > 1
                }
            }
            return dist
        }

        function reverse_complement(seq) {
            # Calculate the reverse complement of a DNA sequence
            n = length(seq)
            rc = ""
            for (i = n; i > 0; i--) {
                base = substr(seq, i, 1)
                if (base == "A") rc = rc "T"
                else if (base == "T") rc = rc "A"
                else if (base == "C") rc = rc "G"
                else if (base == "G") rc = rc "C"
            }
            return rc
        }

        BEGIN {
            # Read ATAC_WHITE into an array
            atac_counter = 0
            while ((getline line < atac_white) > 0) {
                atac_array[++atac_counter] = line
            }
            close(atac_white)

            gex_counter = 0
            # Read GEX_WHITE into an array
            while ((getline line < gex_white) > 0) {
                gex_map[atac_array[++gex_counter]] = line
            }
            close(gex_white)
        }

        NR % 4 == 2 {
            # Extract the key from the sequence
            key = substr($0, length($0) - 15, 16)
            seen[key]++
        }

        END {
            # Sort keys in descending order of counts
            m = asorti(seen, sorted_keys, "@val_num_desc")

            for (j = 1; j <= m; j++) {
                key = sorted_keys[j]
                count = seen[key]
                rc_key = reverse_complement(key)

                # Check for an exact match in preloaded atac_map
                if (rc_key in gex_map) {
                    final_map[rc_key] += count
                    continue
                }

                approx_found = 0
                match_key = ""

                # Check previously processed keys in final_map for approximate matches
                for (prev_key in final_map) {
                    if (final_map[prev_key] > 1 && hamming_distance(rc_key, prev_key) == 1) {
                        approx_found++
                        if (approx_found > 1) break
                        match_key = prev_key
                    }
                }

                if (approx_found > 1) continue

                if (approx_found == 1) {
                    final_map[match_key] += count
                    continue
                }

                # Check for approximate matches in atac_map
                for (atac_key in gex_map) {
                    if (hamming_distance(atac_key, rc_key) == 1) {
                        approx_found++
                        if (approx_found > 1) break
                        match_key = atac_key
                    }
                }

                if (approx_found == 1) {
                    final_map[match_key] += count
                }
            }

            # Print results
            for (key in final_map) {
                print gex_map[key] ", " final_map[key]
            }
        }' | \
        sort -t, -k2,2nr > $OUTPUT_FOLDER/${NAME_PART}/${NAME_PART}_Wei_unfiltered.csv
    }

    export -f process_file_atac
    export OUTPUT_FOLDER ATAC_CB_GREP ATAC_READ2_GREP ATAC_WHITE GEX_WHITE ATAC_FOLDERS

    # Use parallel to process files in parallel, limiting to MAX_PARALLEL_JOBS
    parallel -j $max_samples process_file_atac ::: "${ATAC_FILE_NAMES[@]}"
fi

chmod 777 $OUTPUT_FOLDER
chmod 777 $OUTPUT_FOLDER/*
