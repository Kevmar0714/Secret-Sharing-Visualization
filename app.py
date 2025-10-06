import streamlit as st
import plotly.graph_objects as go
import numpy as np

from logic.shamir import make_shares, lagrange_interpolate_at_zero, gf

st.set_page_config(page_title="Shamir Visualizer", layout="wide")
st.title("Shamir's Secret Sharing (GF(p))")

#--- Controls ----
col1, col2 = st.columns(2)
with col1:
    p = st.number_input("Prime p", min_value=257,value=7919, step=2)
    secret = st.number_input("Secret s (0 ≤ s < p)", min_value=0, value=1234,max_value=p-1)
    n = st.slider("Numer of shares n", 3, 20, 8)
with col2:
    t = st.slider("Threshold t", 2, n, 3)
    reveal_k = st.slider("Shares revealed for reconstruction",1,n,2)
    seed = st.number_input("Random seed", min_value=0,value=0)

np.random.seed(seed)

# --Compute shares ---
GF, f, xs, ys = make_shares(secret=secret, n=n, t=t,p=p)

#pick revealed shares
rx, ry = xs[:reveal_k], ys[:reveal_k]
recon = lagrange_interpolate_at_zero(rx,ry,GF)

st.markdown(f"**Reconstructed f(0):** `{recon}` "f"{'✅ correct' if reveal_k >= t and recon == secret else '⚠️ not guaranteed'}")

# --- Prepare data for plotting (cast GF elements to ints) ---
X_all = [int(x) for x in xs]
Y_all = [int(y) for y in ys]
X_rev = [int(x) for x in rx]
Y_rev = [int(y) for y in ry]
dense_x = GF.Range(1, min(p, 400))
dense_y = f(dense_x)
X_curve = [int(x) for x in dense_x]
Y_curve = [int(y) for y in dense_y]

fig = go.Figure()
fig.add_trace(go.Scatter(x=X_curve, y=Y_curve, mode="lines", name="f(x)"))
fig.add_trace(go.Scatter(x=X_all, y=Y_all, mode="markers", name="All shares"))
fig.add_trace(go.Scatter(x=X_rev, y=Y_rev, mode="markers",marker=dict(size=12, symbol="star"),name="Revealed shares"))
fig.update_layout(height=500, xaxis_title="x (GF element as int)", yaxis_title="y")
st.plotly_chart(fig, use_container_width=True)



