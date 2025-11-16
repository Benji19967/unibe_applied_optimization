"""
Much cleaner plotting using a df and seaborn
"""

from pathlib import Path

import matplotlib.pyplot as plt
import polars as pl
import seaborn as sns

HERE = Path(__file__).resolve().parent
ASSETS_DIR = HERE / "assets"
TIMINGS_DIR = ASSETS_DIR / "timings"

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


def create_df(run_times, run_times_evaluation, num_iterations, obj_func_values):
    g = pl.DataFrame({"grid_size": GRID_SIZES})
    s = pl.DataFrame({"spring_type": SPRING_TYPES})
    m = pl.DataFrame({"method": METHODS})

    df = g.join(s, how="cross").join(m, how="cross")
    df = pl.concat(
        [
            df,
            pl.DataFrame({"run_time": run_times}),
            pl.DataFrame({"run_time_evaluation": run_times_evaluation}),
            pl.DataFrame({"num_iteration": num_iterations}),
            pl.DataFrame({"obj_func_value": obj_func_values}),
        ],
        how="horizontal",
    )

    return df


def plot_objective_function_value(df):
    sns.catplot(
        df,
        x="obj_func_value",
        y="method",
        hue="spring_type",
        col="grid_size",
        kind="bar",
    )
    filename = ASSETS_DIR / "objective_function_value"
    plt.savefig(filename)

    plt.show()


def plot_run_time(df):
    sns.catplot(
        df,
        x="run_time",
        y="method",
        hue="spring_type",
        col="grid_size",
        kind="bar",
    )
    filename = ASSETS_DIR / "run_time"
    plt.savefig(filename)

    plt.show()


def plot_num_iterations(df: pl.DataFrame):
    sns.catplot(
        df,
        x="num_iteration",
        y="method",
        hue="spring_type",
        col="grid_size",
        kind="bar",
    )
    filename = ASSETS_DIR / "num_iterations"
    plt.savefig(filename)

    plt.show()


def plot_percent_eval_of_total(df: pl.DataFrame):
    df = df.with_columns(
        percent_run_time_eval=pl.col("run_time_evaluation") / pl.col("run_time")
    )
    sns.catplot(
        df,
        x="percent_run_time_eval",
        y="method",
        hue="spring_type",
        col="grid_size",
        kind="bar",
    )
    filename = ASSETS_DIR / "percent_eval_of_total"
    plt.savefig(filename)

    plt.show()


def main():
    run_times = []
    run_times_evaluation = []
    obj_func_values = []
    num_iterations = []
    for grid_size in GRID_SIZES:
        for spring_type in SPRING_TYPES:
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

                num_iters = iter[-1]
                obj_func_value = obj[-1]
                run_times.append(total_time)
                run_times_evaluation.append(total_time_evaluation)
                obj_func_values.append(obj_func_value)
                num_iterations.append(num_iters)
    df = create_df(run_times, run_times_evaluation, num_iterations, obj_func_values)
    return df


if __name__ == "__main__":
    df = main()
    print(df)
    # plot_objective_function_value(df)
    # plot_run_time(df)
    # plot_num_iterations(df)
    # plot_percent_eval_of_total(df)
    # df.write_csv(ASSETS_DIR / "summary.csv")
    print(df.describe())
