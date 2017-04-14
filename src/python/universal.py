"""Subtroutines needed by all the systems."""

def get_HNF_diagonals(n):
    """Finds the diagonals of the HNF that reach the target n value.
    
    Args:
        n (int): The target determinant for the HNF.
        
    Retruns:
        diags (list of lists): The allowed values of the determinant.
    """
    
    diags = []
    for i in range(1,n+1):
        if not n%i == 0:
            continue
        else:
            q = n/i
            for j in range(1,q+1):
                if not q%j == 0:
                    continue
                else:
                    diags.append([i,j,q/j])
                    
    return diags
