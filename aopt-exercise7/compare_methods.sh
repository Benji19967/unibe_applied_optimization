get_spring_name() {
    case "$1" in
        2) echo "withLength" ;;
        1) echo "withLength" ;;
        0) echo "noLength" ;;
    esac
}

MAX_ITERATIONS=1000000
ASSETS_DIR=assets

### Gradient Descent
METHOD=gradient_descent

SPRING_SYSTMS_DIR=${ASSETS_DIR}/spring_systems/${METHOD}
TIMINGS_DIR=${ASSETS_DIR}/timings/${METHOD}
mkdir -p ${SPRING_SYSTMS_DIR}
mkdir -p ${TIMINGS_DIR}

for GRID_SIZE in 5 10 20; do
    for SPRING_TYPE in 0 1; do
        SPRING_NAME=$(get_spring_name ${SPRING_TYPE})
        FILENAME=spring_${SPRING_NAME}_${GRID_SIZE}x${GRID_SIZE} 
        ./build/Build/bin/GradientDescent ${SPRING_TYPE} 0 ${GRID_SIZE} ${GRID_SIZE} ${MAX_ITERATIONS} \
            ${SPRING_SYSTMS_DIR}/${FILENAME} \
            &> ${TIMINGS_DIR}/${FILENAME}.txt
    done
done

### Standard Newton
METHOD=standard_newton
NEWTON_METHOD=0 

SPRING_SYSTMS_DIR=${ASSETS_DIR}/spring_systems/${METHOD}
TIMINGS_DIR=${ASSETS_DIR}/timings/${METHOD}
mkdir -p ${SPRING_SYSTMS_DIR}
mkdir -p ${TIMINGS_DIR}

for GRID_SIZE in 5 10 20; do
    SPRING_TYPE=0  # without length
    SPRING_NAME=$(get_spring_name ${SPRING_TYPE})
    FILENAME=spring_${SPRING_NAME}_${GRID_SIZE}x${GRID_SIZE} 
    ./build/Build/bin/NewtonMethods ${NEWTON_METHOD} ${SPRING_TYPE} ${GRID_SIZE} ${GRID_SIZE} ${MAX_ITERATIONS} \
        ${SPRING_SYSTMS_DIR}/${FILENAME} \
        &> ${TIMINGS_DIR}/${FILENAME}.txt
done

for GRID_SIZE in 5 10 20; do
    SPRING_TYPE=2  # with length: use pos. semidefinite H
    SPRING_NAME=$(get_spring_name ${SPRING_TYPE})
    FILENAME=spring_${SPRING_NAME}_${GRID_SIZE}x${GRID_SIZE} 
    ./build/Build/bin/NewtonMethods ${NEWTON_METHOD} ${SPRING_TYPE} ${GRID_SIZE} ${GRID_SIZE} ${MAX_ITERATIONS} \
        ${SPRING_SYSTMS_DIR}/${FILENAME} \
        &> ${TIMINGS_DIR}/${FILENAME}.txt
done

### Projected Newton
METHOD=projected_newton
NEWTON_METHOD=1

SPRING_SYSTMS_DIR=${ASSETS_DIR}/spring_systems/${METHOD}
TIMINGS_DIR=${ASSETS_DIR}/timings/${METHOD}
mkdir -p ${SPRING_SYSTMS_DIR}
mkdir -p ${TIMINGS_DIR}

for GRID_SIZE in 5 10 20; do
    for SPRING_TYPE in 0 1; do
        SPRING_NAME=$(get_spring_name ${SPRING_TYPE})
        FILENAME=spring_${SPRING_NAME}_${GRID_SIZE}x${GRID_SIZE} 
        ./build/Build/bin/NewtonMethods ${NEWTON_METHOD} ${SPRING_TYPE} ${GRID_SIZE} ${GRID_SIZE} ${MAX_ITERATIONS} \
            ${SPRING_SYSTMS_DIR}/${FILENAME} \
            &> ${TIMINGS_DIR}/${FILENAME}.txt
    done
done

### L-BFGS
METHOD=l_bfgs
M=20  # 3-20

SPRING_SYSTMS_DIR=${ASSETS_DIR}/spring_systems/${METHOD}
TIMINGS_DIR=${ASSETS_DIR}/timings/${METHOD}
mkdir -p ${SPRING_SYSTMS_DIR}
mkdir -p ${TIMINGS_DIR}

for GRID_SIZE in 5 10 20; do
    for SPRING_TYPE in 0 1; do
        SPRING_NAME=$(get_spring_name ${SPRING_TYPE})
        FILENAME=spring_${SPRING_NAME}_${GRID_SIZE}x${GRID_SIZE} 
        ./build/Build/bin/LBFGS ${SPRING_TYPE} ${GRID_SIZE} ${GRID_SIZE} ${M} ${MAX_ITERATIONS} \
            ${SPRING_SYSTMS_DIR}/${FILENAME} \
            &> ${TIMINGS_DIR}/${FILENAME}.txt
    done
done

### Gauss-Newton
METHOD=gauss_newton

SPRING_SYSTMS_DIR=${ASSETS_DIR}/spring_systems/${METHOD}
TIMINGS_DIR=${ASSETS_DIR}/timings/${METHOD}
mkdir -p ${SPRING_SYSTMS_DIR}
mkdir -p ${TIMINGS_DIR}

for GRID_SIZE in 5 10 20; do
    for SPRING_TYPE in 0 1; do
        SPRING_NAME=$(get_spring_name ${SPRING_TYPE})
        FILENAME=spring_${SPRING_NAME}_${GRID_SIZE}x${GRID_SIZE} 
        ./build/Build/bin/GaussNewton ${SPRING_TYPE} ${GRID_SIZE} ${GRID_SIZE} ${MAX_ITERATIONS} \
            ${SPRING_SYSTMS_DIR}/${FILENAME} \
            &> ${TIMINGS_DIR}/${FILENAME}.txt
    done
done
