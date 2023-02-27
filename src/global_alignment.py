import sys
from linear_alignment import linear_alignment_int
from affine_alignment import affine_alignment_int
import os

def main():
    
    print("Global alignment: do you want to align two sequences linearly (l) or using an affine gap cost function (a)?")
   
    type_alig = input().strip()
    while type_alig.lower() != 'l' and type_alig.lower() != 'a':
        print("Error. Please introduce a valid input (l/a)")
        type_alig = input().strip()
        
    if type_alig.lower() == 'l':
        linear_alignment_int()
    elif type_alig.lower() == 'a':
        affine_alignment_int()
    


if __name__ == '__main__':
    main()