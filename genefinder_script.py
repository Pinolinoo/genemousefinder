
import matplotlib.pyplot as plt
import seaborn as sns
import functions as ft

genename = "Snca" # type gene name here, first letter capital, rest small

df_small, ages, classes = ft.load_data(genename) #load data from loom file (http://mousebrain.org/downloads.html)

counts_per_age_vector = ft.calculate_total_counts_per_age(df_small, ages, genename)

## make plots of total Counts per Age

fig, axes = plt.subplots(1, 2)

sns.set(font_scale = 1, rc={"figure.figsize":(13, 8)})
sns.lineplot(ax = axes[0], data=counts_per_age_vector, x = "Ages in Embryonic Days", y = "Total Counts of %s" %genename)

sns.set(font_scale = 1, rc={"figure.figsize":(13, 8)})
sns.lineplot(ax = axes[1], data=counts_per_age_vector, x = "Ages in Embryonic Days", y = "%s Counts Normalized by Cells per Age" %genename)

### make heatmap of counts per Class per Age

counts_per_class_per_age  = ft.get_counts_per_class_per_age(df_small, ages, classes, genename)

plt.figure(figsize=(60,20))
heat_map = sns.heatmap(counts_per_class_per_age, linewidth = 0.1 , annot = True, vmax = 1000)
plt.title( "HeatMap for %s Counts per Cluster per Age" %genename )
plt.show()