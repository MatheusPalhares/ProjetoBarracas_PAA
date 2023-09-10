#TeteuCraft e PaulinBaleia

import numpy as np
from scipy.spatial import ConvexHull
import matplotlib.pyplot as plt
import PySimpleGUI as sg

h , w = input().split()
#print(f"Quantidade de homens {h} e quantidade de mulheres {w}")
arr1_a = []
arr1_b = []
arr3_a = []
arr3_b = []
arr5_a = []
arr5_b = []
arr7_a = []
arr7_b = []
#print("Digite os valores de X e Y das barracas masculinas.")
""" def split_input(input, arr):
    coords = input.split()
    for a in coords:
        arr.append(int(a))
    return arr """
for i in range(int(h)):
    inp = []
    inp.clear()
    inp = input().split()
    #print(f"input: {inp}")
    #print("Digite o valor de X." + str(i) + "do canto inferior esquerdo.")
    arr1_a.append(int(inp[0]))
    #print("Digite o valor de Y." + str(i) + "do canto inferior esquerdo.")
    arr1_b.append(int(inp[1]))
    #print("Digite o valor de X." + str(i) + "do canto superior direito.")
    arr3_a.append(int(inp[2]))
    #print("Digite o valor de Y." + str(i) + "do canto superior direito.")
    arr3_b.append(int(inp[3]))
#print("Digite os valores de X e Y das barracas femininas.")
for i in range(int(w)):
    inp2 = []
    inp2.clear()
    inp2 = input().split()
    #print("Digite o valor de X." + str(i) + "do canto inferior esquerdo.")
    arr5_a.append(int(inp2[0]))
    #print("Digite o valor de Y." + str(i) + "do canto inferior esquerdo.")
    arr5_b.append(int(inp2[1]))
    #print("Digite o valor de X." + str(i) + "do canto superior direito.")
    arr7_a.append(int(inp2[2]))
    #print("Digite o valor de Y." + str(i) + "do canto superior direito.")
    arr7_b.append(int(inp2[3]))

aux1 = []
aux2 = []
for i in range(int(h)):
    aux1.append([(arr1_a[i], arr1_b[i]), (arr1_a[i], arr3_b[i]), (arr3_a[i], arr1_b[i]), (arr3_a[i], arr3_b[i])])
    square_edgesM = np.array([[arr1_a[i], arr1_b[i]], [arr1_a[i], arr3_b[i]], [arr3_a[i], arr1_b[i]], [arr3_a[i], arr3_b[i]]])
    hull = ConvexHull(square_edgesM)
    plt.plot(square_edgesM[:,0], square_edgesM[:,1], 'o', color='blue')
    for simplex in hull.simplices:
        plt.plot(square_edgesM[simplex, 0], square_edgesM[simplex, 1], color='blue')
for i in range(int(w)):
    aux2.append([(arr5_a[i], arr5_b[i]), (arr5_a[i], arr7_b[i]), (arr7_a[i], arr5_b[i]), (arr7_a[i], arr7_b[i])])
    square_edgesW = np.array([[arr5_a[i], arr5_b[i]], [arr5_a[i], arr7_b[i]], [arr7_a[i], arr5_b[i]], [arr7_a[i], arr7_b[i]]])
    hull2 = ConvexHull(square_edgesW)
    plt.plot(square_edgesW[:,0], square_edgesW[:,1], 'o', color='red')
    for simplex in hull2.simplices:
        plt.plot(square_edgesW[simplex, 0], square_edgesW[simplex, 1], color='red')
#points = np.random.rand(10, 2)
#print(points)

# Compute the convex hull
#hull4 = ConvexHull(aux2)

np_aux1 = np.array(aux1)
aux1_final = np.concatenate(np_aux1)
np_aux2 = np.array(aux2)
aux2_final = np.concatenate(np_aux2)

hull_M = ConvexHull(aux1_final)
hull_W = ConvexHull(aux2_final)

print(hull_M.vertices)
print(hull_W.vertices)

for simplex in hull_M.simplices:
        plt.plot(aux1_final[simplex, 0], aux1_final[simplex, 1], color='blue')
for simplex in hull_W.simplices:
        plt.plot(aux2_final[simplex, 0], aux2_final[simplex, 1], color='red')

# Check if they are overlapping
""" overlapping = True
for hull in [hull_M, hull_W]:
  # Loop over the edges of the hull
  for i in range(hull.nsimplex):
    # Get the normal vector of the edge
    p1 = hull.points[hull.simplices[i, 0]]
    p2 = hull.points[hull.simplices[i, 1]]
    normal = np.array([p2[1] - p1[1], p1[0] - p2[0]])
    # Project all points of both hulls onto the normal vector
    projection1 = np.dot(hull_M.points, normal)
    projection2 = np.dot(hull_W.points, normal)
    # Check if the projections are disjoint
    if max(projection1) < min(projection2) or max(projection2) < min(projection1):
      # Found a separating line
      overlapping = False
      break
  if not overlapping:
    break """

# Define a small tolerance value
tol = -1

# Check if they are overlapping
overlapping = True
for hull in [hull_M, hull_W]:
  # Loop over the edges of the hull
  for i in range(hull.nsimplex):
    # Get the normal vector of the edge
    p1 = hull.points[hull.simplices[i, 0]]
    p2 = hull.points[hull.simplices[i, 1]]
    normal = np.array([p2[1] - p1[1], p1[0] - p2[0]])
    # Project all points of both hulls onto the normal vector
    projection1 = np.dot(hull_M.points, normal)
    projection2 = np.dot(hull_W.points, normal)
    # Check if the projections are disjoint with some tolerance
    if max(projection1) < min(projection2) - tol or max(projection2) < min(projection1) - tol:
      # Found a separating line
      overlapping = False
      break
  if not overlapping:
    break

#is_separable = not overlapping
plt.show()

layout = [[sg.Text(f"É possível separar? \n\n {not overlapping}")]]

window = sg.Window("display", layout)

while True:
    event, values = window.read()
    if event == sg.WIN_CLOSED:
        break


window.close()

""" root = tk.Tk()
label = tk.Label(root, textvariable=overlapping)
label.pack()
root.mainloop() """

#print(overlapping) # True or False

#print(hull3)
#print(hull4)
# Plot the points and the hull
""" plt.plot(square_edgesW[:,0], square_edgesW[:,1], 'o')
plt.plot(square_edgesM[:,0], square_edgesM[:,1], 'o')
for simplex in hull.simplices:
    plt.plot(square_edgesW[simplex, 0], square_edgesW[simplex, 1], 'k-')
for simplex in hull2.simplices:
    plt.plot(square_edgesM[simplex, 0], square_edgesM[simplex, 1], 'k-')
"""