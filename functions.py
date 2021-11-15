import numpy as np
import loompy
import pandas as pd


def load_data(genename):
    ds_dev = loompy.connect("dev_all.agg.loom", validate=False)

    CountData = pd.DataFrame(data=ds_dev[:, :], index=ds_dev.ra.Gene, columns=ds_dev.ca.ClusterName)

    genename_counts_per_cluster = CountData.loc[genename, :]  # get mean counts per cell from specified Gene per Cluster

    ages = ds_dev.ca.keys()[0:20]  # get ages from Data
    ages = ages[-4:] + ages[:-4]  # right order of Ages

    classes = pd.DataFrame(ds_dev.ca.Class, index=ds_dev.ca.ClusterName)  # get classes of cells

    metadata = {}

    # Make Metadata DataFrame
    for i in ages:
        new = {i: ds_dev.ca[i]}
        metadata.update(new)
    metadata = pd.DataFrame(metadata, index=ds_dev.ca.ClusterName).T
    df = metadata.append(CountData)

    # Make DataFrame with all the Data + Metadata
    df_small = metadata.append(genename_counts_per_cluster)
    df_small = df_small.append(classes.T)
    df_small = df_small.T

    return df_small, ages, classes


def calculate_total_counts_per_age(df_small, ages, genename):
    # total counts
    counts_per_age_vector = np.zeros(len(ages))
    counts_per_cluster_per_age = np.zeros([len(ages), df_small.shape[0]])
    number_of_cells_per_age = np.zeros(len(ages))
    j = 0
    for i in ages:
        for k in range(df_small.shape[0]):
            counts_per_cluster_per_age[j, k] = df_small[i][k] * df_small[genename][k]
            counts_per_age_vector[j] += counts_per_cluster_per_age[j, k]
            number_of_cells_per_age[j] += df_small[i][k]
        j += 1
    ages_number = [7, 8, 8.5, 9, 10, 11, 12, 12.5, 13, 13.5, 14, 14.5, 15, 15.5, 16, 16.25, 16.5, 17, 17.5, 18]
    counts_per_age_vector = pd.DataFrame(counts_per_age_vector, index=ages, columns=["Total Counts of %s" % genename])
    age_df = pd.DataFrame(ages_number, index=ages, columns=["Ages in Embryonic Days"])
    counts_per_age_vector = pd.concat([counts_per_age_vector, age_df], axis=1)
    # normalized per cells per age
    counts_per_age_normalized = np.zeros(len(ages))
    j = 0
    for i in ages:
        counts_per_age_normalized[j] = counts_per_age_vector["Total Counts of %s" % genename][j] / \
                                       number_of_cells_per_age[j]
        j += 1
    counts_per_age_normalized = pd.DataFrame(counts_per_age_normalized, index=ages,
                                             columns=["%s Counts Normalized by Cells per Age" % genename])
    counts_per_age_vector = pd.concat([counts_per_age_vector, counts_per_age_normalized], axis=1)

    return counts_per_age_vector


def get_counts_per_class_per_age(df_small, ages, classes, genename):
    classes_unique = classes[0].unique()
    counts_per_class_per_age = np.zeros([classes_unique.shape[0], len(ages)])
    counts_per_class_per_age = pd.DataFrame(counts_per_class_per_age, index=classes_unique, columns=ages).T

    j = 0
    for i in ages:
        for k in range(df_small.shape[0]):
            current_class = classes[0][k]
            counts_per_class_per_age[current_class][j] += df_small[i][k] * df_small[genename][k]
        j += 1

    return counts_per_class_per_age