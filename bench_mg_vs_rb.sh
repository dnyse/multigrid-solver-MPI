#!/bin/bash -l
#SBATCH --job-name=weak_scale_mg_rb
#SBATCH --output=weak_scale_%j.txt
#SBATCH --partition=singlenode
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=72
#SBATCH --time=06:00:00
#SBATCH --export=NONE
#SBATCH --cpu-freq=2400000-2400000:performance

unset SLURM_EXPORT_ENV

module load likwid intel intelmpi
export I_MPI_PIN=1
export I_MPI_DEBUG=0

# =============================================================================
# WEAK SCALING: MG vs RB-SOR
#
# Weak scaling = constant LOCAL grid per process, increase processes.
#
# Canal:  elongated domain (ratio ~3:1:1), scale primarily in i-direction
#         Original: 200x50x50 with domain 30x4x4
#         Base local grid: 96x32x32 per process (~98k cells)
#         dx=0.3125, dy=dz=0.125 kept constant by scaling domain with grid
#
# Dcavity: cubic domain (ratio 1:1:1), scale uniformly
#          Original: 128x128x128 with domain 1x1x1
#          Base local grid: 48x48x48 per process (~110k cells)
#          dx=dy=dz=0.02083 kept constant
#
# All local dimensions divisible by 16 -> supports up to 4 MG levels.
# Domain lengths scaled proportionally with grid so dx/dy/dz stay constant.
# =============================================================================

WORK_DIR="${SLURM_SUBMIT_DIR}"
cd "${WORK_DIR}" || exit 1

BASE_DIR="${WORK_DIR}/3D-mpi"
EXE_MG="${WORK_DIR}/exe-MG"
EXE_RB="${WORK_DIR}/exe-RB"
FILENAME="${WORK_DIR}/weak_scaling_results.csv"

# --- Weak scaling configurations ---
# Format: "ranks imax jmax kmax xlength ylength zlength"
#
# CANAL: elongated, scale i-direction first (long axis of pipe)
# Base: 96x32x32 cells, domain 30x4x4
# dx=30/96=0.3125, dy=4/32=0.125, dz=4/32=0.125
#
# Scaling strategy: grow the pipe longer first, then widen
#   ranks  decomp     global grid      domain
#   1      1x1x1      96x32x32         30x4x4
#   2      2x1x1      192x32x32        60x4x4
#   4      4x1x1      384x32x32        120x4x4
#   8      4x2x1      384x64x32        120x8x4
#   16     4x2x2      384x64x64        120x8x8
#   32     8x2x2      768x64x64        240x8x8
#   64     8x4x2      768x128x64       240x16x8

CANAL_CONFIGS=(
    "1   96  32  32   30.0  4.0  4.0"
    "2   192 32  32   60.0  4.0  4.0"
    "4   384 32  32   120.0 4.0  4.0"
    "8   384 64  32   120.0 8.0  4.0"
    "16  384 64  64   120.0 8.0  8.0"
    "32  768 64  64   240.0 8.0  8.0"
    "64  768 128 64   240.0 16.0 8.0"
)

# DCAVITY: cubic, scale uniformly
# Base: 48x48x48 cells, domain 1x1x1
# dx=dy=dz=1/48=0.02083
#
# Scaling strategy: grow all directions as evenly as possible
#   ranks  decomp     global grid      domain
#   1      1x1x1      48x48x48         1x1x1
#   2      2x1x1      96x48x48         2x1x1
#   4      2x2x1      96x96x48         2x2x1
#   8      2x2x2      96x96x96         2x2x2
#   16     4x2x2      192x96x96        4x2x2
#   32     4x4x2      192x192x96       4x4x2
#   64     4x4x4      192x192x192      4x4x4

DCAVITY_CONFIGS=(
    "1   48  48  48   1.0  1.0  1.0"
    "2   96  48  48   2.0  1.0  1.0"
    "4   96  96  48   2.0  2.0  1.0"
    "8   96  96  96   2.0  2.0  2.0"
    "16  192 96  96   4.0  2.0  2.0"
    "32  192 192 96   4.0  4.0  2.0"
    "64  192 192 192  4.0  4.0  4.0"
)

# Solver configurations
MG_LEVELS=(2 3)

# Short simulation time — enough for stable averages, not too long
BENCH_TE="5.0"

# Repeats for statistics
REPEATS=3

# --- Functions ---

setup_config() {
    local config_file=$1
    local imax=$2
    local jmax=$3
    local kmax=$4
    local xlength=$5
    local ylength=$6
    local zlength=$7
    local mg_levels=$8

    sed -i "s/^imax[[:space:]]\+.*$/imax          ${imax}/" "${BASE_DIR}/${config_file}"
    sed -i "s/^jmax[[:space:]]\+.*$/jmax          ${jmax}/" "${BASE_DIR}/${config_file}"
    sed -i "s/^kmax[[:space:]]\+.*$/kmax          ${kmax}/" "${BASE_DIR}/${config_file}"
    sed -i "s/^xlength[[:space:]]\+.*$/xlength       ${xlength}/" "${BASE_DIR}/${config_file}"
    sed -i "s/^ylength[[:space:]]\+.*$/ylength       ${ylength}/" "${BASE_DIR}/${config_file}"
    sed -i "s/^zlength[[:space:]]\+.*$/zlength       ${zlength}/" "${BASE_DIR}/${config_file}"
    sed -i "s/^te[[:space:]]\+.*$/te       ${BENCH_TE}/" "${BASE_DIR}/${config_file}"
    sed -i "s/^levels[[:space:]]\+.*$/levels        ${mg_levels}/" "${BASE_DIR}/${config_file}"
}

run_single() {
    local ranks=$1
    local solver_type=$2
    local config_file=$3
    local imax=$4
    local jmax=$5
    local kmax=$6
    local xlength=$7
    local ylength=$8
    local zlength=$9
    local mg_levels=${10}
    local run_id=${11}

    local executable
    if [ "${solver_type}" == "MG" ]; then
        executable="${EXE_MG}"
    else
        executable="${EXE_RB}"
    fi

    local problem_name
    problem_name=$(basename "${config_file}" .par)

    # Setup config
    setup_config "${config_file}" "${imax}" "${jmax}" "${kmax}" \
                 "${xlength}" "${ylength}" "${zlength}" "${mg_levels}"

    # Pin processes
    local np_1=$(($ranks - 1))
    export I_MPI_PIN_PROCESSOR_LIST=0-$np_1

    # Run with timeout (10 min per run to catch hangs)
    local OUTPUT
    OUTPUT=$(timeout 600 mpirun -n ${ranks} "${executable}" "${BASE_DIR}/${config_file}" 2>&1)
    local exit_code=$?

    # Parse results
    local WALLTIME="FAIL"
    local ITERATIONS="FAIL"
    local RESIDUAL="FAIL"

    if [ ${exit_code} -eq 124 ]; then
        WALLTIME="TIMEOUT"
        ITERATIONS="TIMEOUT"
        RESIDUAL="TIMEOUT"
    elif [ ${exit_code} -ne 0 ]; then
        WALLTIME="CRASH"
        ITERATIONS="CRASH"
        RESIDUAL="CRASH"
    else
        local result
        result=$(echo "${OUTPUT}" | grep "Solution took")
        if [ -n "$result" ]; then
            WALLTIME=$(echo "${result}" | grep -oP 'Solution took \K[0-9.]+')
            ITERATIONS=$(echo "${result}" | grep -oP 'in \K[0-9]+')
        fi

        local res_line
        res_line=$(echo "${OUTPUT}" | grep "Final residiuum")
        if [ -n "$res_line" ]; then
            RESIDUAL=$(echo "${res_line}" | grep -oP '[0-9]+\.[0-9]+')
        fi
    fi

    # Cells per process
    local total_cells=$(( imax * jmax * kmax ))
    local cells_per_proc=$(( total_cells / ranks ))

    local lvl_str
    if [ "${solver_type}" == "MG" ]; then
        lvl_str="${mg_levels}"
    else
        lvl_str="N/A"
    fi

    # CSV
    echo "${problem_name},${solver_type},${lvl_str},${ranks},${imax}x${jmax}x${kmax},${xlength}x${ylength}x${zlength},${cells_per_proc},${ITERATIONS},${WALLTIME},${RESIDUAL},${run_id}" >> "${FILENAME}"

    # Console
    printf "  %-7s %-3s lvl=%-3s ranks=%-3s grid=%-15s -> %6s iters %9ss res=%s (run %d)\n" \
        "${problem_name}" "${solver_type}" "${lvl_str}" "${ranks}" \
        "${imax}x${jmax}x${kmax}" \
        "${ITERATIONS}" "${WALLTIME}" "${RESIDUAL}" "${run_id}"
}

# =============================================================================
# MAIN
# =============================================================================

if [ ! -d "${BASE_DIR}" ]; then
    echo "ERROR: Directory ${BASE_DIR} does not exist!"
    exit 1
fi

# Backup configs
for prob in canal dcavity; do
    if [ -f "${BASE_DIR}/${prob}.par" ]; then
        cp "${BASE_DIR}/${prob}.par" "${BASE_DIR}/${prob}.par.backup"
    fi
done

# --- Compile ---
echo "=========================================="
echo "COMPILATION"
echo "=========================================="

make -C "${BASE_DIR}" distclean 2>/dev/null || true

echo "--- Building MG solver ---"
make -C "${BASE_DIR}" SOLVER=mg TAG=ICX
if [ $? -ne 0 ]; then echo "ERROR: MG compile failed"; exit 1; fi
cp "${BASE_DIR}/exe-ICX" "${EXE_MG}"

make -C "${BASE_DIR}" distclean
echo "--- Building RB solver ---"
make -C "${BASE_DIR}" SOLVER=rb TAG=ICX
if [ $? -ne 0 ]; then echo "ERROR: RB compile failed"; exit 1; fi
cp "${BASE_DIR}/exe-ICX" "${EXE_RB}"

# --- Setup results ---
rm -f "${FILENAME}"
echo "Problem,Solver,MGLevels,Ranks,Grid,Domain,CellsPerProc,Iterations,Walltime_s,FinalResidual,RunID" > "${FILENAME}"

# --- Print config ---
echo ""
echo "=========================================="
echo "WEAK SCALING BENCHMARK"
echo "=========================================="
echo "Sim time:  ${BENCH_TE}s"
echo "MG levels: ${MG_LEVELS[*]}"
echo "Repeats:   ${REPEATS}"
echo ""
echo "Canal (elongated ~3:1:1, local ~96x32x32 = 98304 cells/proc):"
for config_str in "${CANAL_CONFIGS[@]}"; do
    read -r r i j k xl yl zl <<< "${config_str}"
    printf "  ranks=%-3s grid=%-15s domain=%-17s local=%-15s\n" \
        "${r}" "${i}x${j}x${k}" "${xl}x${yl}x${zl}" \
        "$((i/r > 0 ? i : i))x${j}x${k}"
done
echo ""
echo "Dcavity (cubic ~1:1:1, local ~48x48x48 = 110592 cells/proc):"
for config_str in "${DCAVITY_CONFIGS[@]}"; do
    read -r r i j k xl yl zl <<< "${config_str}"
    printf "  ranks=%-3s grid=%-15s domain=%-17s\n" \
        "${r}" "${i}x${j}x${k}" "${xl}x${yl}x${zl}"
done
echo "=========================================="
echo ""

# --- Run Canal ---
echo "=== CANAL (elongated domain, Re=100) ==="
echo ""

for config_str in "${CANAL_CONFIGS[@]}"; do
    read -r ranks imax jmax kmax xlength ylength zlength <<< "${config_str}"
    local_cells=$(( imax * jmax * kmax / ranks ))

    echo "--- Ranks=${ranks}, Grid=${imax}x${jmax}x${kmax}, Domain=${xlength}x${ylength}x${zlength} (${local_cells} cells/proc) ---"

    # RB-SOR
    for run in $(seq 1 ${REPEATS}); do
        run_single "${ranks}" "RB" "canal.par" "${imax}" "${jmax}" "${kmax}" \
                   "${xlength}" "${ylength}" "${zlength}" "1" "${run}"
    done

    # MG with different levels
    for lvl in "${MG_LEVELS[@]}"; do
        for run in $(seq 1 ${REPEATS}); do
            run_single "${ranks}" "MG" "canal.par" "${imax}" "${jmax}" "${kmax}" \
                       "${xlength}" "${ylength}" "${zlength}" "${lvl}" "${run}"
        done
    done

    echo ""
done

# --- Run Dcavity ---
echo "=== DCAVITY (cubic domain, Re=1000) ==="
echo ""

for config_str in "${DCAVITY_CONFIGS[@]}"; do
    read -r ranks imax jmax kmax xlength ylength zlength <<< "${config_str}"
    local_cells=$(( imax * jmax * kmax / ranks ))

    echo "--- Ranks=${ranks}, Grid=${imax}x${jmax}x${kmax}, Domain=${xlength}x${ylength}x${zlength} (${local_cells} cells/proc) ---"

    # RB-SOR
    for run in $(seq 1 ${REPEATS}); do
        run_single "${ranks}" "RB" "dcavity.par" "${imax}" "${jmax}" "${kmax}" \
                   "${xlength}" "${ylength}" "${zlength}" "1" "${run}"
    done

    # MG with different levels
    for lvl in "${MG_LEVELS[@]}"; do
        for run in $(seq 1 ${REPEATS}); do
            run_single "${ranks}" "MG" "dcavity.par" "${imax}" "${jmax}" "${kmax}" \
                       "${xlength}" "${ylength}" "${zlength}" "${lvl}" "${run}"
        done
    done

    echo ""
done

# --- Restore configs ---
for prob in canal dcavity; do
    if [ -f "${BASE_DIR}/${prob}.par.backup" ]; then
        mv "${BASE_DIR}/${prob}.par.backup" "${BASE_DIR}/${prob}.par"
    fi
done

# --- Summary ---
echo "=========================================="
echo "BENCHMARK COMPLETE"
echo "=========================================="
echo "Results saved to: ${FILENAME}"
echo ""
echo "Results:"
column -t -s',' "${FILENAME}"
