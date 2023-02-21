#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
import random
import math
import matplotlib.pyplot as plt
import matplotlib.colors as colors

############## Base functions

def initialization_with_iridophores(length, width):
    '''
    Initialize a matrix of dimension length*length with zeros and 
    an horizontal "band" of iridophores in the matrix center
    Parameters
    ----------
    length : int
        Number of lines/columns of the matrix
    width : int
        The width of the iridophores band.
        
    Returns
    -------
    array
        A matrix which represents the zebrafish lattice
    int
        First line of iridophores band
    int
        Last line of iridophores band
    -------
    '''
    mat = np.zeros((length,length), dtype=int)
    irid_start = (length//2) - (width//2)
    irid_end = irid_start + width
    return(mat,irid_start,irid_end)    
           

def random_neighbour_of(row,col,length,mat):
    '''
    Returns a random neighbor of the cell mat[row,col]
    Parameters
    ----------
    row : int
        Currrent cell row
    col : int
        Current cell column
    length : int
        Number of lines and columns of the matrix.
    mat : array
        Matrix which represents the lattice.

    Returns
    -------
    int
        The value of one of the 8 neighbors
    '''    
    vertical = random.choice([(row-1)%length,(row+1)%length,row])
    #if the neighbour row is the same as the cell, we cannot chose the same column as the cell
    if vertical == row: 
        horizontal = random.choice([(col-1)%length,(col+1)%length])
    else: 
        horizontal = random.choice([(col-1)%length,(col+1)%length,col])
    return (mat[vertical,horizontal])

def display(mat,irid_start,irid_end,h,lx):
    '''
  Prints the matrix with colors instead of numbers.
  Empty cell : 0 -> white
  Xantophore : 1 -> yellow
  Melanophore : 2 -> black
  '''
    plt.figure()
    plt.pcolormesh(mat,cmap=colors.ListedColormap(['white','yellow','black']))
    if irid_end != 0:
        plt.axhline(y=irid_start,color='red')
        plt.axhline(y=irid_end,color='red')
    plt.title(f" h = {h} and lx = {lx}")
    plt.show()
    

############## Event functions: each function is related to one of the seven events
  

def birth_xantophore(cell) :
    '''   
    First, checks if the cell is empty. Then, if it is the case, a xantophore will be born in this cell.
    ----------
    cell : int
        The value of the current cell

    Returns
    -------
    int
        1 if a xantophore is born/if the cell was already occupied by a xantophore
        2 if the cell was already occupied by a melanophore
    '''
    if cell == 0:
        return 1
    return cell

def birth_melanophore(row,col,cell,irid_start,irid_end) :
    '''  
    First, checks if the cell is empty, and is not located on top of the iridophores band. 
    Then, if it respects the conditions, a melanophore will be born in this cell.
    Parameters
    ----------
    row : int
        Currrent cell row
    col : int
        Current cell column
    cell : int
        The value of the current cell
    irid_start : int
        First line of iridophores band.
    irid_end : int
        Last line of iridophores band.

    Returns
    -------
    int
        2 if a melanophore is born/if the cell was already occupied by a melanophore
        1 if the cell was already occupied by a xantophore
        0 if the cell is on top of the iridophores band
    '''
    if cell == 0 and (row < irid_start or row > irid_end):
        return 2
    return cell

def death_xantophore(cell):
    '''
    If the cell is a xantophore, it dies
    ----------
    cell : int
        The value of the current cell

    Returns
    -------
    int
        0 if a xantophore is dead
        2 if the cell was occupied by a melanophore
    '''
    if cell == 1:
        return 0
    return cell 

def death_melanophore(cell):
    '''
    If the cell is a melanophore, it dies
    ----------
    cell : int
        The value of the current cell

    Returns
    -------
    int
        0 if a melanophore is dead
        1 if the cell was occupied by a xantophore

    '''
    if cell == 2:
        return 0
    return cell


def short_range_competitivity_death_xanto(row,col,mat,cell,length) :
    '''
    Fifth event of the algorithm, its represents melanophores inhibition effect on xantophores.
    If the cell is a xantophore, and one the random neighbour is a melanophore, the xantophore dies.
    For borders of the matrix, a neighbour can be the cell on the opposite border.
     Parameters
    ----------
    row : int
        Currrent cell row
    col : int
        Current cell column
    mat : TYPE
        Matrix which represents the lattice.
    cell : int
        The value of the current cell
    length : int
        Number of lines and columns of the matrix.

    Returns
    -------
    int
        0 if a xantophore is dead
        1 if the neighbour of the xantophore was not a melanophore
        2 if the cell was occupied by a melanophore
    '''
    if cell == 1:
        neighbour = random_neighbour_of(row,col,length,mat)
        if neighbour == 2:
            return 0
    return cell
        
 

def short_range_competitivity_death_melano(row,col,mat,cell,length) :
    '''
    Fifth event of the algorithm, its represents xantophores inhibition effect on melanophores.
    If the cell is a melanophore, and one the random neighbour is a xantophore, the melanophore dies.
    For borders of the matrix, a neighbour can be the cell on the opposite border.


     Parameters
    ----------
    row : int
        Currrent cell row
    col : int
        Current cell column
    mat : array
        Matrix which represents the lattice.
    cell : int
        The value of the current cell
    length : int
        Number of lines and columns of the matrix.

    Returns
    -------
    int
        0 if a melanophore is dead
        1 if the neighbour of the melanophore was not a xantophore
        2 if the cell was occupied by a xantophore

    '''
    if cell == 2:
        neighbour = random_neighbour_of(row,col,length,mat)
        if neighbour == 1:
            return 0
    return cell

def long_range_effect_birth_melano(row,col,mat,cell,h,length,irid_start,irid_end) :
    '''
    An random neighbour at a distance h from the current cell is selected. 
    If this neighbour is a xantophore and the current cell is empty, 
    a melanophore will be born on this cell

     Parameters
    ----------
    row : int
        Currrent cell row
    col : int
        Current cell column
    mat : array
        Matrix which represents the lattice.
    cell : int
        The value of the current cell
    h : int
        The distance between the xantophore and its neighbour
    length : int
        Number of lines and columns of the matrix.
    irid_start : int
        First line of iridophores band.
    irid_end : int
        Last line of iridophores band.

    Returns
    -------
        int
        2 if a melanophore is born/if the cell was already occupied by a melanophore
        1 if the cell was already occupied by a xantophore
        0 if the cell is the selected neighbour was not a xantophore

    '''
    # obtenir l'angle teta entre 0 et 2pi
    teta = random.uniform(0, (2 * math.pi))
    # les coordonnées du voisin à distance h (on en choisit qu'un seul)
    delta_x = round(h * math.cos(teta))
    delta_y = round(h * math.sin(teta))
    if cell == 0 and (row < irid_start or row > irid_end):
        neighbour = mat[(row+delta_y)%length,(col+delta_x)%length]
        if neighbour == 1:
            return 2
    return cell

def run(lx,h,irid):
    '''
    Initializes a matrix which represents zebrafish lattice, eventually with an iridophores band
    Simulates the process according to parameters. Multiple iterations (default = 10⁹) are made.
    For each one, a random cell is selected, 
    and a random event occures (according to probabilities given in the article).
    At the end, the matrix is converted to a plot, which represents the lattice with colors.

    Parameters
    ----------
    lx : int
        Strength of enhancement effect of xantophores on melanophores
    h : int
        The required distance for melanophore to enhance xantophores birth
    int
        The width of the iridophores band.
    

    '''
    bx = 1 # birth rate
    dx = 0 # death rate
    sx = 1 # xantophores strength
    
    bm = 0 # birth rate
    dm = 0 # death rate
    sm = 1 # melanophores strength
    
    length = 100 # number of lines/columns of the matrix
    nb_iter = 1000000000 #number of iterations
    
    # Matrix initialization
    if irid == 0:
        mat = np.zeros((length,length),dtype=int)
        #if there is not iridophores band, we define irid_start and irid_end, as the boundaries of the matrix
        irid_start = length
        irid_end = 0
    
    else:    
        (mat,irid_start,irid_end) = initialization_with_iridophores(length, irid)
    
    # Events rates
    sum_rates = bx + bm + dx + dm + sm + sx + lx
    pBX = bx/sum_rates
    pBM = bm/sum_rates
    pDX = dx/sum_rates
    pDM = dm/sum_rates
    pCDX = sm/sum_rates 
    pCDM = sx/sum_rates
    
    for time in range (nb_iter) :
        row = random.randint(0, length-1)
        col = random.randint(0, length-1)
        cell = mat[row,col]
        Pk = random.random() #probabilité de l'événement pour cette itération
        
           
        if Pk <= pBM : 
            mat[row,col] = birth_melanophore(row,col,cell,irid_start,irid_end)
            
        elif Pk <= pBX + pBM : 
            mat[row,col] = birth_xantophore(cell)
        
        elif Pk <= pBX + pBM + pDX : 
            mat[row,col] = death_xantophore(cell)
            
        elif Pk <= pBX + pBM + pDX + pDM : 
            mat[row,col] = death_melanophore(cell)
            
        elif Pk <= pBX + pBM + pDX + pDM + pCDX : 
            mat[row,col] = short_range_competitivity_death_xanto(row,col,mat,cell,length)
            
        elif Pk <= pBX + pBM + pDX + pDM + pCDX + pCDM : 
            mat[row,col] = short_range_competitivity_death_melano(row,col,mat,cell,length)
            
        else: 
            mat[row,col] = long_range_effect_birth_melano(row,col,mat,cell, h, length,irid_start,irid_end)
            
        if time%100000000 == 0:            
            print(f"{time} iterations")
            
    display(mat,irid_start,irid_end,h,lx)        
    
#Function calling    
run(lx=2.5,h=16,irid=1)
