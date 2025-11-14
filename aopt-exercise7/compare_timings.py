import csv
from pathlib import Path

import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent
TIMINGS_DIR = HERE / "assets" / "timings"

METHODS = [
    "gauss_newton",
    "gradient_descent",
    "l_bfgs",
    "projected_newton",
    "standard_newton",
]
GRID_SIZES = [
    5,
    10,
    20,
]
SPRING_TYPES = [
    "noLength",
    "withLength",
]


def parse_timings(filename: Path):
    total_time = None
    total_time_evaluation = None
    iter = []
    obj = []
    stopping_criterion = []

    with open(filename, "r") as f:
        for line in f:
            if line.startswith("iter"):
                line_list = line.split()
                iter.append(int(line_list[1]))
                obj.append(float(line_list[4]))
                stopping_criterion.append(float(line_list[7]))
            if line.startswith("total time  "):
                line_list = line.split()
                total_time = float(line_list[3][:-1])
            if line.startswith("total time evaluation"):
                line_list = line.split()
                total_time_evaluation = float(line_list[4][:-1])

    return total_time, total_time_evaluation, iter, obj, stopping_criterion


def plot_objective_function_value(iter, obj, method, grid_size, spring_type):
    _, ax = plt.subplots()
    ax.set_title(
        f"{method}: obj function value vs num iterations {grid_size}x{grid_size} {spring_type}"
    )
    ax.plot(iter, obj)
    filename = (
        HERE
        / "assets"
        / "obj_function_values"
        / f"{method}_spring_{spring_type}_{grid_size}x{grid_size}"
    )
    plt.savefig(filename)


def plot_run_time(methods, total_times, grid_size, spring_type):
    # Function to add centered value labels
    def add_labels(x, y):
        for i in range(len(x)):
            plt.text(i, y[i], y[i], ha="center")  # Aligning text at center

    _, ax = plt.subplots()
    ax.bar(methods, total_times, color=["#00979D"])

    add_labels(methods, total_times)
    plt.xticks(rotation=10)
    ax.set_title(f"Total run time per method {grid_size}x{grid_size} {spring_type}")
    ax.set_ylabel("Total run time (s)")
    filename = (
        HERE
        / "assets"
        / "total_run_times"
        / f"spring_{spring_type}_{grid_size}x{grid_size}"
    )
    plt.savefig(filename)


def plot_num_iterations(methods, num_iterations, grid_size, spring_type):
    _, ax = plt.subplots()
    ax.bar(methods, num_iterations, color=["#00979D"])
    plt.xticks(rotation=10)

    ax.set_title(f"Num iterations per method {grid_size}x{grid_size} {spring_type}")
    filename = (
        HERE
        / "assets"
        / "num_iterations"
        / f"spring_{spring_type}_{grid_size}x{grid_size}"
    )
    plt.savefig(filename)


def main():
    summary_filename = HERE / "assets" / "summary.csv"
    with open(summary_filename, "w+") as f:
        csv.writer(f).writerow(
            (
                "grid_side_size",
                "spring_element_type",
                "method",
                "num_iterations",
                "total_time",
                "total_time_evaluation",
                "obj_func_value",
            )
        )
        for grid_size in GRID_SIZES:
            for spring_type in SPRING_TYPES:
                run_times = []
                obj_func_values = []
                num_iterations = []
                for method in METHODS:
                    filename = (
                        TIMINGS_DIR
                        / method
                        / f"spring_{spring_type}_{grid_size}x{grid_size}.txt"
                    )
                    (
                        total_time,
                        total_time_evaluation,
                        iter,
                        obj,
                        stopping_criterion,
                    ) = parse_timings(filename)
                    plot_objective_function_value(
                        iter, obj, method, grid_size, spring_type
                    )

                    num_iters = iter[-1]
                    obj_func_value = obj[-1]
                    csv.writer(f).writerow(
                        (
                            grid_size,
                            spring_type,
                            method,
                            num_iters,
                            total_time,
                            total_time_evaluation,
                            obj_func_value,
                        )
                    )
                    run_times.append(total_time)
                    obj_func_values.append(obj_func_value)
                    num_iterations.append(num_iters)
                plot_run_time(METHODS, run_times, grid_size, spring_type)
                plot_num_iterations(METHODS, num_iterations, grid_size, spring_type)


if __name__ == "__main__":
    main()
