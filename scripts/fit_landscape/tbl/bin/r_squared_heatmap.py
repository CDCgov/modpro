import pandas as pd
from itertools import combinations
import statsmodels.api as sm
import dash
from dash import dcc, html
import plotly.express as px
from sklearn.linear_model import LinearRegression

df = pd.read_csv("scw_bvic_01-22-2024.tsv", sep="\t", header=0)
variants = df["variant_hash"]  # storing the variants to use for metadata later
columns_to_drop = [  # remove any non-numeric or unwanted columns
    "subtype",
    "protein",
    "clade",
    "variant_hash",
    "substitution_list",
    "surveillance_period",
    "aa_seq",
    "data_source",
]
df = df.drop(columns=columns_to_drop)


def linear_regression_sklearn(data, x, y):
    model = LinearRegression().fit(data[x].values.reshape(-1, 1), data[y])
    r_squared = model.score(data[x].values.reshape(-1, 1), data[y])
    return model, r_squared


# Get all pairs of unique columns excluding "variant_hash"
column_pairs = [(x, y) for x, y in combinations(df.columns, 2)]

# Create a DataFrame to store R-squared values
heatmap_data = pd.DataFrame(index=df.columns, columns=df.columns)

# Populate the DataFrame with R-squared values
for pair in column_pairs:
    x, y = pair
    model, r_squared = linear_regression_sklearn(df, x, y)
    heatmap_data.loc[x, y] = r_squared
    heatmap_data.loc[y, x] = r_squared

# Convert R-squared values to numeric
heatmap_data = heatmap_data.apply(pd.to_numeric)

df["variant_hash"] = variants

# Initialize Dash app
app = dash.Dash(__name__)

# Create layout
app.layout = html.Div(
    [
        dcc.Graph(
            id="heatmap",
            figure=px.imshow(
                heatmap_data, x=heatmap_data.columns, y=heatmap_data.index
            ),
        ),
        dcc.Graph(id="scatter-plot"),
    ]
)


# Callback to update scatter plot based on selected cell in the heatmap
@app.callback(
    dash.dependencies.Output("scatter-plot", "figure"),
    [dash.dependencies.Input("heatmap", "clickData")],
)
def display_scatter_plot(clickData):
    if clickData is None:
        # If no cell is clicked, return an empty figure
        return px.scatter()

    # Get selected columns from clickData
    x, y = clickData["points"][0]["x"], clickData["points"][0]["y"]

    # Perform linear regression
    model, _ = linear_regression_sklearn(df, x, y)

    # Create scatter plot with the line of best fit
    scatter_fig = px.scatter(
        df, x=x, y=y, title=f"Scatter Plot: {x} vs {y}", hover_data=["variant_hash"]
    )
    scatter_fig.add_traces(
        px.line(x=df[x], y=model.predict(df[x].values.reshape(-1, 1))).data
    )

    return scatter_fig


# Run the app
if __name__ == "__main__":
    app.run_server(debug=True)
