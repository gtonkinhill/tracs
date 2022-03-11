import plotly.graph_objects as go
import plotly.offline as offline
import numpy as np
from collections import Counter
from datetime import datetime


def plot_heatmap(dates, clusters, outfile):

    cluster_date_count = Counter()
    for i in dates:
        cluster_date_count[(dates[i][0], clusters[i])] += 1

    cluster_counts = []
    dates = []
    cluster = []
    for dc in cluster_date_count:
        cluster_counts.append(cluster_date_count[dc])
        dates.append(datetime.fromisoformat(dc[0]))
        cluster.append(str(dc[1] + 1))

    fig = go.Figure(data=go.Heatmap(
        z=cluster_counts, x=dates, y=cluster, colorscale='Viridis'))

    fig.update_layout(title='Genomes per transmission cluster per day',
                      yaxis_nticks=len(set(cluster)) + 1)

    fig.update_xaxes(title_text='Time')
    fig.update_yaxes(title_text='Transmission Cluster')

    offline.plot(fig, filename=outfile, auto_open=True)

    return