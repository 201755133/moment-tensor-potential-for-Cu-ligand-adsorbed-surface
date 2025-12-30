import matplotlib.pyplot as plt
import numpy as np

fig = plt.figure()
ax = fig.add_subplot()

E =  [0.00, -0.74, -2.78, -3.25, -4.22, -4.34, -5.54, -4.60] 
E1 = [0.00, -2.35, -3.23, -4.42, -5.55, -6.84, -7.02, -6.17]

#ax.set_ylim([min(E+E1)-0.2,max(E+E1)+0.2])

ax.plot([0.1,1.1],[E[0],E[0]],'-',color = 'grey',linewidth = 4)
ax.plot([1.1,1.4],[E[0],E[1]],':',color = 'grey')
ax.plot([1.4,2.4],[E[1],E[1]],'-',color = 'grey',linewidth = 4)
ax.plot([2.4,2.7],[E[1],E[2]],':',color = 'grey')
ax.plot([2.7,3.7],[E[2],E[2]],'-',color = 'grey',linewidth = 4)
ax.plot([3.7,4.0],[E[2],E[3]],':',color = 'grey')
ax.plot([4.0,5.0],[E[3],E[3]],'-',color = 'grey',linewidth = 4)
ax.plot([5.0,5.3],[E[3],E[4]],':',color = 'grey')
ax.plot([5.3,6.3],[E[4],E[4]],'-',color = 'grey',linewidth = 4)
ax.plot([6.3,6.6],[E[4],E[5]],':',color = 'grey')
ax.plot([6.6,7.7],[E[5],E[5]],'-',color = 'grey',linewidth = 4)
# Energy of CH2CHO
ax.plot([7.7,8.0],[E[5],E[6]],':',color = 'grey')
ax.plot([8.0,9.0],[E[6],E[6]],'-',color = 'grey',linewidth = 4)
# Energy of O*
ax.plot([7.7,8.0],[E[5],E[7]],':',color = 'grey')
ax.plot([8.0,9.0],[E[7],E[7]],'-',color = 'grey',linewidth = 4)
# last step is ethanol
#ax.plot([9.0,9.3],[E[6],E[8]],':',color = 'grey')
#ax.plot([9.3,10.3],[E[8],E[8]],'-',color = 'grey',linewidth = 4)

# last step is liquied water
#ax.plot([9.0,9.3],[E[7],E[9]],':',color = 'grey')
#ax.plot([9.3,10.3],[E[9],E[9]],'-',color = 'grey',linewidth = 0)


ax.plot([0.1,1.1],[E1[0],E1[0]],'-',color = 'blue',linewidth = 4)
ax.plot([1.1,1.4],[E1[0],E1[1]],':',color = 'blue')
ax.plot([1.4,2.4],[E1[1],E1[1]],'-',color = 'blue',linewidth = 4)
ax.plot([2.4,2.7],[E1[1],E1[2]],':',color = 'blue')
ax.plot([2.7,3.7],[E1[2],E1[2]],'-',color = 'blue',linewidth = 4)
ax.plot([3.7,4.0],[E1[2],E1[3]],':',color = 'blue')
ax.plot([4.0,5.0],[E1[3],E1[3]],'-',color = 'blue',linewidth = 4)
ax.plot([5.0,5.3],[E1[3],E1[4]],':',color = 'blue')
ax.plot([5.3,6.3],[E1[4],E1[4]],'-',color = 'blue',linewidth = 4)
ax.plot([6.3,6.6],[E1[4],E1[5]],':',color = 'blue')
ax.plot([6.6,7.7],[E1[5],E1[5]],'-',color = 'blue',linewidth = 4)
# Energy of CH2CHO
ax.plot([7.7,8.0],[E1[5],E1[6]],':',color = 'blue')
ax.plot([8.0,9.0],[E1[6],E1[6]],'-',color = 'blue',linewidth = 4)
# Energy of O*
ax.plot([7.7,8.0],[E1[5],E1[7]],':',color = 'blue')
ax.plot([8.0,9.0],[E1[7],E1[7]],'-',color = 'blue',linewidth = 4)
# last step is ethanol
#ax.plot([9.0,9.3],[E1[6],E1[8]],':',color = 'blue')
#ax.plot([9.3,10.3],[E1[8],E1[8]],'-',color = 'blue',linewidth = 4)
# last step is liquied water
#ax.plot([9.0,9.3],[E1[7],E1[9]],':',color = 'blue')
#ax.plot([9.3,10.3],[E1[9],E1[9]],'-',color = 'blue',linewidth = 0)


#ax.annotate('0V, pH 0 Electrode',xy = [0.1,min(E+E1) - 0.1],color = 'grey',  fontsize=12)
#ax.annotate('-0.59V, pH 0 Electrode', xy = [0.1,min(E+E1) - 0.8], color = 'blue',  fontsize=12)

#an = ax.annotate('CO$(gas)',xy = [0.1,E[0]+0.04], color = 'grey')
#an.draggable()
#an = ax.annotate('*COCOH',xy = [1.4,E[1]+0.04], color = 'grey')
#an.draggable()
#an = ax.annotate('*C$_2$O',xy = [2.7,E[2]+0.04], color = 'grey')
#an.draggable()
#an = ax.annotate('*CHCO',xy = [4.0,E[3]+0.04], color = 'grey')
#an.draggable()
#an = ax.annotate('*CHCHO',xy = [5.3,E[4]+0.04], color = 'grey')
#an.draggable()
#an = ax.annotate('*CH$_2$CHO',xy = [6.6,E[5]+0.04], color = 'grey')
#an.draggable()
#an = ax.annotate('*CH$_3$CHO',xy = [8.0,E[6]+0.04], color = 'grey')
#an.draggable()
#an = ax.annotate('*O + C$_2$H$_4$(gas)',xy = [8.0,E[7]+0.04], color = 'grey')
#an.draggable()
#an = ax.annotate('*CH$_3$CH$_2$OH',xy = [9.3,E[8]+0.04], color = 'grey')
#an.draggable()
#an = ax.annotate('*H$_2$O',xy = [9.3,E[9]+0.04], color = 'grey')
#an.draggable()

# At applied potential
#an = ax.annotate('CO$(gas)',xy = [0.1,E1[0]+0.04], color = 'blue')
#an.draggable()
#an = ax.annotate('*COCOH',xy = [1.4,E1[1]+0.04], color = 'blue')
#an.draggable()
#an = ax.annotate('*C$_2$O',xy = [2.7,E1[2]+0.04], color = 'blue')
#an.draggable()
#an = ax.annotate('*CHCO',xy = [4.0,E1[3]+0.04], color = 'blue')
#an.draggable()
#an = ax.annotate('*CHCHO',xy = [5.3,E1[4]+0.04], color = 'blue')
#an.draggable()
#an = ax.annotate('*CH$_2$CHO',xy = [6.6,E1[5]+0.04], color = 'blue')
#an.draggable()
#an = ax.annotate('*CH$_3$CHO',xy = [8.0,E1[6]+0.04], color = 'blue')
#an.draggable()
#an = ax.annotate('*O + C$_2$H$_4$(gas)',xy = [8.0,E1[7]+0.04], color = 'blue')
#an.draggable()
#an = ax.annotate('*CH$_3$CH$_2$OH',xy = [9.3,E1[8]+0.04], color = 'blue')
#an.draggable()
#an = ax.annotate('*H$_2$O',xy = [9.3,E1[9]+0.04], color = 'blue')
#an.draggable()

# Set y-axis limits
ax.set_ylim([-8, 6])
#ax.set_xlim([0, 12])
ax.set_xlim([0, 11])
ax.set_xticks([])
plt.xlabel("Reaction Coordinate", fontsize = 15)
plt.ylabel("$\Delta G$(eV)", fontsize = 15)

# Increase the size of tick labels
ax.tick_params(axis='both', which='major', labelsize=15)  # For major ticks
ax.tick_params(axis='both', which='minor', labelsize=12)  # For minor ticks (if needed)

# Save the plot as a PDF file
plt.savefig("free_energy_diagram.png", format='png')
plt.show()


