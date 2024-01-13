import igraph
import matplotlib.pyplot as plt

# Paramètres du modèle binomial
S0 = 100  # Prix initial de l'actif sous-jacent
u = 1.1   # Facteur d'augmentation du prix à chaque étape
d = 0.9   # Facteur de diminution du prix à chaque étape
n_steps = 2  # Nombre d'étapes

# Création d'un objet graphique (arbre)
g = igraph.Graph.Tree(n_steps + 1,int(n_steps * (u + d) / 2) , mode=1)

# Renommer les nœuds pour refléter les prix
for i in range(len(g.vs)):
    g.vs[i]["name"] = f"Prix: {S0 * (u ** i) * (d ** (n_steps - i)):.2f}"

# Visualisation de l'arbre
layout = g.layout_reingold_tilford(root=[0])
igraph.plot(g, layout=layout, bbox=(400, 300), margin=20)

plt.scatter([0,1,2,3],[5,6,7,8])