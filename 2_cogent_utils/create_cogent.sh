#!/bin/sh

usage() {
  cat <<EOF
Usage: $(basename "$0") arg1

Creates a cogent generator script

Available options:

arg1            Input file that contains the fasta files to be transformed
arg2            Name of the generated script's output folder
-h, --help      Print this help and exit

EOF
  exit
}


parse_params() {
  # default values of variables set from params
  flag=0
  input_file="${1-GenomeAll}"
  output_folder="${2-CoGent_ID}"

  while :; do
    case "${1-}" in
    -h | --help)
      usage
      ;;
    -?*)
      die "Unknown option: $1"
      ;;
    *)
      break
      ;;
    esac
    shift
  done

  args=("$@")

  return 0
}


parse_params "$@"            

printf "%s\n\n" \#\!/bin/sh

printf "%s\n\n"  "mkdir ${output_folder}"

cat "${input_file}" | while read line; do
	echo $line > temp
	printf "%s" "awk '/>/{printf(\"+"
	printf "%s" "`awk '{print $2}' temp;`"
	printf "%s" "-01-%006d \", i++)}1'  "
	printf "%s" "`awk '{print $1}' temp;`"
	printf "%s" " > ${output_folder}/";
	printf "%s" "`awk '{print $2}' temp;`"
	printf "%s \n" ".faa;"

	printf "%s" "sed -i 's/>//g' ${output_folder}/";
	printf "%s" "`awk '{print $2}' temp;`";
	printf "%s \n" ".faa;";

	printf "%s" "sed -i 's/+/>/g' ${output_folder}/";
	printf "%s" "`awk '{print $2}' temp;`";
	printf "%s \n" ".faa;";
done;
rm temp
