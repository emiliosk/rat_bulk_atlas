import pandas as pd
import plotly.express as px


def make_sunburst_plot(df, color_col, color_key_col, path_cols, title=None, save_path=None):
  """
    Create an interactive sunburst plot from a dataframe.
    
    Args:
        df (pd.DataFrame): Input dataframe
        color_col (str): Column containing color values
        color_key_col (str): Column to use as key for color mapping
        path_cols (list): List of columns defining hierarchy from outer to inner
        title (str, optional): Plot title
        save_path (str, optional): Path to save HTML output
        
    Returns:
        plotly.graph_objects.Figure: The sunburst plot figure
    """
  # Sort dataframe by path columns
  df = df.sort_values(path_cols).reset_index(drop=True)
  
  # Create color mapping
  color_map = df[[color_key_col, color_col]].drop_duplicates().set_index(color_key_col).to_dict()[color_col]
  
  # Create sunburst plot
  fig = px.sunburst(df, 
                    path=path_cols,
                    title=title,
                    color=color_key_col,
                    color_discrete_map=color_map)
  
  # Save if path provided
  if save_path:
    fig.write_html(save_path)

  return fig


sample_meta = pd.read_csv("./local/integrated_data/filtered_samples_obs.tsv", sep = '\t')

fig = make_sunburst_plot(
  df=sample_meta,
  color_col="organ_color",
  color_key_col="organ_name",
  path_cols=["organ_name", "consensus_tissue_name", "tissue_name", "sample_id"],
  title="Pig: metadata",
  save_path="./figures/sunburst_plot_w_sample.html"
)

fig.write_image("./figures/sunburst_plot_w_sample.pdf")

fig = make_sunburst_plot(
  df=sample_meta,
  color_col="organ_color",
  color_key_col="organ_name",
  path_cols=["organ_name", "consensus_tissue_name", "tissue_name"],
  title="Pig: metadata",
  save_path="./figures/sunburst_plot.html"
)

fig.write_image("./figures/sunburst_plot.pdf")

# make_sunburst_plot(
#     df=sample_meta,
#     color_col="organ_color",
#     color_key_col="organ_name",
#     path_cols=["organ_name", "consensus_tissue_name", "region_tissue_name", "tissue_name"],
#     title="Mouse:metadata",
#     save_path="./figures/sunburst_plot_w_region.html"
# )
