import matplotlib
import hail as hl

matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn


def histogram(ht:hl.Table, location:str, plot_path:str=None, plot_tmp_path:str= '/home/danfengc/tmp.png', bins:int=100):
    """
        This function is to draw histogram of a column of hail.Table through converting it to
        pandas dataframe first, and then use matplotlib module.
        Note
        ----
        :param Table ht: Input Hail Table
        :param str location: column storing data for drawing histogram
        :param str plot_path: google bucket for saving the output figure
        :param str plot_tmp_path: path for saving the temporary figure in the master node
        :param int bins: number of bins to be plotted
        """
    df = ht.to_pandas()
    print(df[location].head())
    hist = df.hist(column=location, bins=bins)
    plt.savefig(plot_tmp_path, dpi=300)
    plt.close()
    if plot_path is not None:
        hl.utils.hadoop_copy('file://%s'%plot_tmp_path, plot_path)


def scatter(ht: hl.Table, x_location:str, y_location:str, color_location:str=None, plot_path:str=None, plot_tmp_path:str= '/home/danfengc/tmp.png'):
    """
        This function is to draw scatter plot using two columns of hail.Table through converting it to
        pandas dataframe first, and then use matplotlib module.
        Note
        ----
        :param Table ht: Input Hail Table
        :param str x_location: column storing data for x-axis of the scatterplot
        :param str y_location: column storing data for y-axis of the scatterplot
        :param str color_location: column storing data for coloring the scatterplot
        :param str plot_path: google bucket for saving the output figure
        :param str plot_tmp_path: path for saving the temporary figure in the master node
        """
    df = ht.to_pandas()
    fg = seaborn.FacetGrid(data=df, hue=color_location, aspect=1, height=6)
    fg.map(plt.scatter, x_location, y_location).add_legend()
    plt.savefig(plot_tmp_path, dpi=300)
    plt.close()
    if plot_path is not None:
        hl.utils.hadoop_copy('file://%s' % plot_tmp_path, plot_path)