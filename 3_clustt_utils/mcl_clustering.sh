#!/usr/bin/env bash

set -Eeuo pipefail

cd "$(dirname "${BASH_SOURCE[0]}")" >/dev/null 2>&1

trap cleanup SIGINT SIGTERM ERR EXIT

usage() {
  cat <<EOF
Usage: $(basename "$0") arg1 arg2 arg3 arg4

Pairs to mcl families pipeline script

Available options:

arg1            Input pairs file
arg2            Name of the mcl families output file
arg3            I mcl input parameter
arg4            Path of the mcl executable
-h, --help      Print this help and exit
-v, --verbose   Print script debug info

EOF
  exit
}

cleanup() {
  trap - SIGINT SIGTERM ERR EXIT
  # script cleanup here
}

setup_colors() {
  if [[ -t 2 ]] && [[ -z "${NO_COLOR-}" ]] && [[ "${TERM-}" != "dumb" ]]; then
    NOCOLOR='\033[0m' RED='\033[0;31m' GREEN='\033[0;32m' ORANGE='\033[0;33m' BLUE='\033[0;34m' PURPLE='\033[0;35m' CYAN='\033[0;36m' YELLOW='\033[1;33m'
  else
    NOCOLOR='' RED='' GREEN='' ORANGE='' BLUE='' PURPLE='' CYAN='' YELLOW=''
  fi
}

msg() {
  echo >&2 -e "${1-}"
}

die() {
  local msg=$1
  local code=${2-1} # default exit status 1
  msg "$msg"
  exit "$code"
}

parse_params() {
  # default values of variables set from params
  flag=0
  input_file="${1-target.pairs}"
  families_output_file="${2-families.output.txt}"
  i_mcl_param="${3-1.8}"
  mcl_executable="${4-/home/ouzounis/bcpl/bin/mcl}"
  graph_file="${input_file%.*}.graph"
  mcl_output_file="out.${graph_file}.I${i_mcl_param}"
  
  while :; do
    case "${1-}" in
    -h | --help)
      usage
      ;;
    -v | --verbose)
      set -x
      ;;
    --no-color)
      NO_COLOR=1
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

  # check required params and arguments
  #[[ ${#args[@]} -eq 0 ]] && die "Missing script arguments"

  return 0
}

CHECK_MARK="\033[0;32m\xE2\x9C\x94\033[0m"

show_spinner()
{
    local -r pid="${1}"
    local -r delay='0.75'
    local spinstr='\|/-'
    local temp
    while ps a | awk '{print $1}' | grep -q "${pid}"; do
        temp="${spinstr#?}"
        printf " [%c]  " "${spinstr}"
        spinstr=${temp}${spinstr%"${temp}"}
        sleep "${delay}"
        #printf "\b\b\b\b\b\b"
    done
    echo -e "\\r${CHECK_MARK} done"
}

parse_params "$@"

echo "${input_file}"

printf "Counting ${input_file}\n"
awk '{print $1}' "${input_file}" | sort | uniq | wc &
show_spinner "$!"
printf "\n\n"

printf "Creating the graph file: %s\n" "${graph_file}"
awk '$(NF-1)<0.005' "${input_file}" | awk '{print $1,$2,$3}' | sort -rn -k 3 > "${graph_file}" &
show_spinner "$!"
printf "\n\n"

printf "Send graph file to mcl\n"
"${mcl_executable}" "${graph_file}" --abc -I "${i_mcl_param}" -o "${mcl_output_file}"
echo -e "\\r${CHECK_MARK} done"
printf "\n\n"

printf "Creating families file: %s from the mcl output file: %s\n" "${families_output_file}" "${mcl_output_file}"
awk '{printf "%s_%s\t%s\n",NR,NF,$0}' "${mcl_output_file}" > "${families_output_file}" &
show_spinner "$!"

printf "\n"

