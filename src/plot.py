from itertools import accumulate
from math import log2
from matplotlib.lines import Line2D
import matplotlib.pyplot as pyplot
import os

def a_lg_b_list (a_list, b_list):
    return [(a * log2(b)) for (a, b) in zip(a_list, b_list)]

def accumulated_MB_list (list_):
    return list(accumulate([(value / 1024 / 1024 / 8) for value in list_]))

def plot_trend_lines(output_base_name, state, alphabet_size_list, grammar_rules_size_list, grammar_trie_size_list, text_size_list):
    trend_line_list = []
    legend_list = []
    color_list = []

    text_size_0_lg_256_in_MB = (text_size_list[0] * log2(256) / 1024 / 1024 / 8)
    text_size_0_lg_alphabet_size_0_in_MB = (text_size_list[0] * log2(alphabet_size_list[0]) / 1024 / 1024 / 8)
    trend_line_size = len(text_size_list)

    trend_line_list.append(accumulated_MB_list(alphabet_size_list))
    trend_line_list.append(accumulated_MB_list(grammar_rules_size_list))
    trend_line_list.append(accumulated_MB_list(grammar_trie_size_list))
    trend_line_list.append(accumulated_MB_list(text_size_list))
    legend_list.extend([r"$\sigma_{d}$", r"$r_{d}$", r"$t_{d}$", r"$n_{d}$"])

    size_lists = [alphabet_size_list[:-1], alphabet_size_list[1:], grammar_rules_size_list[1:], grammar_trie_size_list[1:], text_size_list[:-1], text_size_list[1:]]
    label_list = [r"$\sigma_{d}$", r"$\sigma_{d + 1}$", r"$r_{d + 1}$", r"$t_{d + 1}$", r"$n_{d}$", r"$n_{d + 1}$"]
    for (a_size_list, a_label) in zip(size_lists, label_list):
        for (b_size_list, b_label) in zip(size_lists, label_list):
            trend_line = accumulated_MB_list(a_lg_b_list(a_size_list, b_size_list))
            trend_line_list.append(trend_line)
            legend_list.append(a_label + " lg " + b_label)

    filtered_trend_line_list = []
    filtered_legend_list = []

    for (index, trend_line) in enumerate(trend_line_list):
            if (trend_line[-1] > text_size_0_lg_256_in_MB):
                if (state == "all"):
                    filtered_trend_line_list.append(trend_line)
                    filtered_legend_list.append(legend_list[index])
                    color_list.append("red")
            elif (trend_line[-1] > text_size_0_lg_alphabet_size_0_in_MB):
                if (state == "all" or state == "warning"):
                    filtered_trend_line_list.append(trend_line)
                    filtered_legend_list.append(legend_list[index])
                    color_list.append("orange")
            else:
                filtered_trend_line_list.append(trend_line)
                filtered_legend_list.append(legend_list[index])
                color_list.append("green")

    (figure_, axes_) = pyplot.subplots(nrows = 1, ncols = 1, figsize = (16, 12))

    for (index, trend_line) in enumerate(filtered_trend_line_list):
        axes_.plot(trend_line, color = color_list[index], linestyle = ":", marker = "${:02d}$".format(index), markersize = 16)

    if (state == "all"):
        axes_.plot([text_size_0_lg_256_in_MB] * trend_line_size, color = "red")
        filtered_legend_list.append(r"$n_{0}$ lg $256$")
    axes_.plot([text_size_0_lg_alphabet_size_0_in_MB] * trend_line_size, color = "orange")
    filtered_legend_list.append(r"$n_{0}$ lg $\sigma_{0}$")

    axes_.set_title("Trend line over depth", fontsize = 20)
    axes_.set_xlabel("depth (d)", fontsize = 20)
    axes_.set_ylabel("accumulated space (MB)", fontsize = 20)
    axes_.set_xticks(list(range(int(trend_line_size * 1.8))))
    axes_.ticklabel_format(axis = 'y', scilimits = (0, 0))
    axes_.tick_params(axis = "both", which = "both", labelsize = 18)
    axes_.legend(filtered_legend_list, ncol = 2, fontsize = 16)
    figure_.savefig(output_base_name + ".trendline." + state + ".png")
    pyplot.close(figure_)

def parse_csv_file (input_root_directory, output_root_directory, file):
    with open(os.path.join(input_root_directory, file), 'r') as input_file:
        lines = input_file.readlines()
        alphabet_size_list = [int(value) for value in lines[0].split(',')]
        grammar_rules_size_list = [int(value) for value in lines[1].split(',')]
        grammar_trie_size_list = [int(value) for value in lines[2].split(',')]
        text_size_list = [int(value) for value in lines[3].split(',')]

    output_base_name = os.path.join(output_root_directory, os.path.splitext(file)[0])
    for state in ["all", "warning", "safe"]:
        plot_trend_lines(output_base_name, state, alphabet_size_list, grammar_rules_size_list, grammar_trie_size_list, text_size_list)

def traverse_csv_files ():
    input_root_directory = "../output"
    output_root_directory = "../output"
    # input_root_directory = "../output/csv"
    # output_root_directory = "../output/figure"
    with os.scandir(input_root_directory) as entries:
        for entry in entries:
            if entry.is_file() and os.path.splitext(entry.name)[1] == ".csv":
                parse_csv_file(input_root_directory, output_root_directory, entry.name)

if __name__ == "__main__":
    traverse_csv_files()
