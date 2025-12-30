import matplotlib.pyplot as plt
import numpy as np

fig = plt.figure()
ax = fig.add_subplot()

E =  [0.00, -0.44,  -0.74, -3.77, -4.99, -4.41] # cuni without CO
E1 = [0.00, -1.47,  -2.35, -5.55, -6.52, -5.97] # cuni with CO


ax.set_ylim([min(E+E1)-0.2,max(E+E1)+0.2])

ax.plot([0.1,1.1],[E[0],E[0]],'-',color = 'grey',linewidth = 4)
ax.plot([1.1,1.4],[E[0],E[1]],':',color = 'grey')

ax.plot([1.4,2.4],[E[1],E[1]],'-',color = 'grey',linewidth = 4)
ax.plot([2.4,2.7],[E[1],E[2]],':',color = 'grey')

ax.plot([2.7,3.7],[E[2],E[2]],'-',color = 'grey',linewidth = 4)
ax.plot([3.7,4.0],[E[2],E[3]],':',color = 'grey')

ax.plot([4.0,5.0],[E[3],E[3]],'-',color = 'grey',linewidth = 4)
ax.plot([5.0,5.3],[E[3],E[4]],':',color = 'grey')
ax.plot([5.3,6.3],[E[4],E[4]],'-',color = 'grey',linewidth = 4)

ax.plot([5.0,5.3],[E[3],E[5]],':',color = 'grey')
ax.plot([5.3,6.3],[E[5],E[5]],'-',color = 'grey',linewidth = 4)




ax.plot([0.1,1.1],[E1[0],E1[0]],'-',color = 'blue',linewidth = 4)
ax.plot([1.1,1.4],[E1[0],E1[1]],':',color = 'blue')

ax.plot([1.4,2.4],[E1[1],E1[1]],'-',color = 'blue',linewidth = 4)
ax.plot([2.4,2.7],[E1[1],E1[2]],':',color = 'blue')

ax.plot([2.7,3.7],[E1[2],E1[2]],'-',color = 'blue',linewidth = 4)
ax.plot([3.7,4.0],[E1[2],E1[3]],':',color = 'blue')

ax.plot([4.0,5.0],[E1[3],E1[3]],'-',color = 'blue',linewidth = 4)
ax.plot([5.0,5.3],[E1[3],E1[4]],':',color = 'blue')
ax.plot([5.3,6.3],[E1[4],E1[4]],'-',color = 'blue',linewidth = 4)


ax.plot([5.0,5.3],[E1[3],E1[5]],':',color = 'blue')
ax.plot([5.3,6.3],[E1[5],E1[5]],'-',color = 'blue',linewidth = 4)


#ax.annotate('CuNi ', xy=(0.1, -0.5), color='grey', fontsize=12)
#ax.annotate('CuNi-*CO', xy=(0.1, -1.0), color='blue', fontsize=12)
#ax.annotate('CuNi ',xy = [0.1,min(E+E1) - 0.1],color = 'grey',  fontsize=12)
#ax.annotate('CuNi-*CO', xy = [0.1,min(E+E1) - 0.8], color = 'blue',  fontsize=12)

#an = ax.annotate('CO',xy = [0.1,E[0]+0.04], color = 'grey')
#an.draggable()
#an = ax.annotate('*COCO*',xy = [1.4,E[1]+0.04], color = 'grey')
#an.draggable()
#an = ax.annotate('*COCOH',xy = [2.7,E[2]+0.04], color = 'grey')
#an.draggable()
#an = ax.annotate('*CH$_2$CHO',xy = [4.0,E[3]+0.04], color = 'grey')
#an.draggable()
#an = ax.annotate('*CCH',xy = [5.3,E[4]+0.04], color = 'grey')
#an.draggable()
#an = ax.annotate('*HOCHCH',xy = [5.3,E[5]+0.04], color = 'grey')
#an.draggable()

#At applied potential
#an = ax.annotate('CO',xy = [0.1,E1[0]+0.04], color = 'blue')
#an.draggable()
#an = ax.annotate('*COCO',xy = [1.4,E1[1]+0.04], color = 'blue')
#an.draggable()
#an = ax.annotate('*COCOH',xy = [2.7,E1[2]+0.04], color = 'blue')
#an.draggable()
#an = ax.annotate('*CH$_2$CHO',xy = [4.0,E1[3]+0.04], color = 'blue')
#an.draggable()
#an = ax.annotate('*CCH',xy = [5.3,E1[4]+0.04], color = 'blue')
#an.draggable()
#an = ax.annotate('*HOCHCH',xy = [5.3,E1[5]+0.04], color = 'blue')
#an.draggable()


# Set y-axis limits
ax.set_ylim([-8, 2])
ax.set_xlim([0, 9])
ax.set_xticks([])
plt.xlabel("Reaction Coordinate", fontsize = 15)
plt.ylabel("$\Delta G$(eV)", fontsize = 15)

# Increase the size of tick labels
ax.tick_params(axis='both', which='major', labelsize=15)  # For major ticks
ax.tick_params(axis='both', which='minor', labelsize=12)  # For minor ticks (if needed)

# Save the plot as a PDF file
plt.savefig("free_energy_diagram.png", format='png')
plt.show()


