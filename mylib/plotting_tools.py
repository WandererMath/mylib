import matplotlib.pyplot as plt


def stacked_hist(labels, values, output):
    # Example input
    # data = {
    #     'Type A': [1, 2, 2, 3, 3, 3],
    #     'Type B': [2, 3, 4, 4],
    #     'Type C': [1, 4, 5]
    # }

    # Extract values and labels
    # values = list(data.values())
    # labels = list(data.keys())

    # Plot stacked histogram
    plt.hist(values, bins=36, stacked=True, label=labels)

    plt.xlabel('Read Length')
    plt.ylabel('Frequency')
    # plt.title('Stacked Histogram by Type')
    plt.legend()
    plt.tight_layout()
    plt.savefig(output)
