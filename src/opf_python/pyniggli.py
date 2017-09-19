"""Contains the main class and functions for the niggli reduced
cell. As described in the papers found at:

https://www.mendeley.com/viewer/?fileId=74bc20b9-a7a5-8e3d-4608-6ce09ea453e0&documentId=7b2a0aec-6dcf-3475-8834-96ddbb760220

https://www.mendeley.com/viewer/?fileId=39464089-ceb6-6b27-6e45-ade432b467cc&documentId=7dbd825c-282f-3f5a-99a3-cc89a1a5a11e

Author: Wiley S. Morgan 2017
"""

import numpy as np

class reduced_cell(object):
    """This class contains the methods necessary to reduce a lattice to
    it's niggli reduced cell.

    Attributes:
        original (numpy ndarray): The original cell vectors.
        niggli (numpy ndarray): The niggli reduced cell vectors.
        C (numpy ndarray): The transformation matrix to transform from the original
            to the niggli cell (O' = O C).
        volume (float): The volume of the cell.        
    
    Examples:
        The following examples show how to generate a niggli reduced cell.

        >>> import numpy an np
        >>> from pyniggli import reduced_cell 
        >>> A = np.transpose([[0.5,0,0.5],[0,3,0],[0.5,0,-0.5]])
        >>> B = reduced_cell(A)
        >>> # The niggli cell can be retrieved with
        >>> print(np.transpose(B.niggli))
            array([[-0.5  0.  -0.5],[-0.5  0.   0.5],[ 0.  -3.   0. ]])
        >>> # To get the transformation matrix for the niggli cell.
        >>> print(np.transpose(B.C))
            array([[-0.5,0,-0.5],[-0.5,0,0.5],[0,-3,0]])
    """

    def __init__(self,A,eps = None):
        """Initial setup of cell.
        Args:
            A (numpy ndarray): A 3 by 3 matrix containing the lattice vectors as columns.
            eps (optional float): Floating point tollerance for comparisons, 
                default is 1E-5.

        Rasise:
            ValueError: if the input is not a 3 by 3 matrix.
            ValueError: if the input has a determinant of 0.
            RuntimeError: if the niggli cell is not found within 100 iterations.
        """

        if not isinstance(A,np.ndarray):
            A = np.array(A)
        if A.shape != (3,3):
            raise ValueError("The input basis must be a 3 by 3 matrix with "
                             "the lattice vectors as columns of the matrix.")
        
        if np.linalg.det(A) ==0:
            raise ValueError("The cell specified has a volume of zero.")
        else:
            self.volume = np.linalg.det(A)
            
        self.original = np.array(A)
        if eps is None:
            self.eps = 1E-5
        else:
            self.eps = eps

        self._niggli_reduction()
        self.niggli = np.dot(self.original,self.C)


    def _niggli_reduction(self):
        """Performs the niggli reduction of the given lattice.

        Raises:
            RuntimeError: if the niggli cell is not found within 100 iterations.
        """

        count, reduced, self.C = 0, False, np.array([[1,0,0],[0,1,0],[0,0,1]]) 

        A = np.dot(self.original[:,0],self.original[:,0])
        B = np.dot(self.original[:,1],self.original[:,1])
        C = np.dot(self.original[:,2],self.original[:,2])
        xi = 2.0 * np.dot(self.original[:,1],self.original[:,2])
        eta = 2.0 * np.dot(self.original[:,2],self.original[:,0])
        zeta = 2.0 *np.dot(self.original[:,0],self.original[:,1])

        while not reduced and count <=100:            
            reduced = True
            count += 1
            #1
            if (A-self.eps)>B or (abs(abs(A)-abs(B))<self.eps and (abs(xi)-self.eps)>abs(eta)):
                A, B = self._swap(A,B)
                xi, eta = self._swap(xi,eta)
                self.C = np.dot(self.C,[[0,-1,0],[-1,0,0],[0,0,-1]])
                
            #2
            if (B-self.eps)>C or (abs(abs(C)-abs(B))<self.eps and (abs(eta)-self.eps)>abs(zeta)):
                B, C = self._swap(B,C)
                eta, zeta = self._swap(eta,zeta)
                self.C = np.dot(self.C,[[-1,0,0],[0,0,-1],[0,-1,0]])
                reduced = False
                continue
                #go to 1
            #3
            if (eta*xi*zeta-self.eps) > 0:
                M = self._find_C3(xi,eta,zeta)
                xi, eta, zeta = abs(xi), abs(eta), abs(zeta)
                self.C = np.dot(self.C,M)
                
            #4
            if not 0 < (eta*xi*zeta-self.eps) and not(xi<-self.eps and eta<-self.eps and zeta<-self.eps):
                M = self._find_C4(xi,eta,zeta)
                xi, eta, zeta = -abs(xi), -abs(eta), -abs(zeta)
                self.C = np.dot(self.C,M)
                
            #5
            if (abs(xi)-self.eps)>B or (np.allclose(xi,B) and 2*eta<(zeta-self.eps)) or (np.allclose(xi,-B) and zeta<(0-self.eps)):
                self.C = np.dot(self.C,np.array([[1,0,0],[0,1,-np.sign(xi)],[0,0,1]]))
                C = B+C-xi*np.sign(xi)
                eta = eta-zeta*np.sign(xi)
                xi = xi-2*B*np.sign(xi)
                reduced = False
                continue
                #go to 1
            #6
            if (abs(eta)-self.eps)>A or (np.allclose(eta,A) and (2*xi<(zeta-self.eps))) or (np.allclose(eta,-A) and zeta<(0-self.eps)):
                self.C = np.dot(self.C,np.array([[1,0,-np.sign(eta)],[0,1,0],[0,0,1]]))
                C = A+C-eta*np.sign(eta)
                xi = xi-zeta*np.sign(eta)
                eta = eta-2*A*np.sign(eta)
                reduced = False
                continue
                #go to 1
            #7
            if (abs(zeta)-self.eps)>A or (np.allclose(zeta,A) and (2*xi<(eta-self.eps))) or (np.allclose(zeta,-A) and eta<(0-self.eps)):
                self.C = np.dot(self.C,np.array([[1,-np.sign(zeta),0],[0,1,0],[0,0,1]]))
                B = A+B-zeta*np.sign(zeta)
                xi = xi-eta*np.sign(zeta)
                zeta = zeta-2*A*np.sign(zeta)
                reduced = False
                continue
                #go to 1
            #8
            if xi+eta+zeta+A+B<(0-self.eps) or (abs(xi+eta+zeta+A+B)<self.eps and (2*(A+eta)+zeta-self.eps)>0):
                C = A+B+C+xi+eta+zeta
                xi = 2*B+xi+zeta
                eta = 2*A+eta+zeta
                self.C = np.dot(self.C,np.array([[1,0,1],[0,1,1],[0,0,1]]))
                reduced = False
                continue
                #go to 1

        if count >= 100:
            raise RuntimeError("Could not reduce the cell in 100 iterations. This could be "
                               "because of floating point error, try providing a smaller eps "
                               "value.")

        if not self._niggli_check(A,B,C,xi,eta,zeta,self.eps): # pragma: no cover
            raise RuntimeError("Cell reduction incorroct A: {0}, B: {1}, C: {2}, "
                               "xi: {3}, eta: {4}, zeta: {5}. Please submit a bug "
                               "report to: https://github.com/wsmorgan/pyniggli/"
                               "issues".format(A,B,C,xi,eta,zeta))


    @staticmethod
    def _niggli_check(A,B,C,xi,eta,zeta,eps):
        """Checks that the niggli reduced cell satisfies the niggli conditions.
        Conditions listed at: https://arxiv.org/pdf/1203.5146.pdf.

        Args:
            A (float): a.a
            B (float): b.b 
            C (float): c.c
            xi (float): 2*b.c
            eta (float): 2*c.a
            zeta (float): 2*a.b
        
        Returns:
            False if niggli conditons aren't met.
        """

        if not (A-eps > 0 and (A < B-eps or np.allclose(A,B,atol=eps)) and
                (B < C-eps or np.allclose(B,C,atol=eps))):
            print("1")
            return False
        
        if np.allclose(A,B,atol=eps) and not (abs(xi) < abs(eta)-eps or
                                              np.allclose(abs(xi),abs(eta),atol=eps)):
            print("2")
            return False
            
        if np.allclose(B,C,atol=eps) and not (abs(eta) < abs(zeta)-eps
                                              or np.allclose(abs(eta),abs(zeta),atol=eps)):
            print("3")
            return False

        if not ((xi-eps > 0 and eta-eps > 0 and zeta-eps > 0) or
                ((xi < 0-eps or np.allclose(xi,0,atol=eps))
                 and (eta < 0-eps or np.allclose(eta,0,atol=eps))
                 and (zeta < 0-eps or np.allclose(zeta,0,atol=eps)))):
            print("4")
            return False

        if not (abs(xi) < B-eps or np.allclose(abs(xi),B,atol=eps)):
            print("5")
            return False 
            
        if not ((abs(eta) < A-eps or np.allclose(abs(eta),A,atol=eps)) and (abs(zeta) < A-eps or
                                                               np.allclose(abs(zeta),A,atol=eps))):
            print("6")
            return False

        if not (C < A+B+C+xi+eta+zeta-eps or np.allclose(C, A+B+C+xi+eta+zeta,atol=eps)):
            print("7")
            return False

        if np.allclose(xi,B,atol=eps) and not (zeta < 2.*eta-eps or
                                               np.allclose(zeta,2.*eta,atol=eps)):
            print("8")
            return False

        if np.allclose(eta,A,atol=eps) and not (zeta < 2.*xi-eps or
                                                np.allclose(zeta,2.*xi,atol=eps)):
            print("9")
            return False

        if np.allclose(zeta,A,atol=eps) and not (eta < 2.*xi-eps or
                                                 np.allclose(eta,2.*xi,atol=eps)):
            print("10")
            return False

        if np.allclose(xi,-B,atol=eps) and not np.allclose(zeta,0,atol=eps):
            print("11")
            return False

        if np.allclose(eta,-A,atol=eps) and not np.allclose(zeta,0,atol=eps):
            print("12")
            return False

        if np.allclose(zeta,-A,atol=eps) and not np.allclose(eta,0,atol=eps):
            print("13")
            return False

        if np.allclose(C,A+B+C+xi+eta+zeta,rtol=0.0) and not ((2.*A+2.*eta+zeta) < 0-eps or
                                                     np.allclose(2.*A+2.*eta+zeta,0,atol=eps)):
            print("14")
            return False

        return True

    @staticmethod
    def _swap(A,B):
        """Swaps the values of A and B.

        Args:
            A (float): A value.
            B (float): Another value.

        Returns:
            B, A (float,float): The values of A and B swapped.
        """

        return B,A

    def _find_C3(self,xi,eta,zeta):
        """Finds the correct transformation matrix given the values of xi, eta, 
        and zeta for step 3.

        Args:
            xi (float): The value of xi.
            eta (float): The value of eta.
            zeta (float): The value of zeta.

        Returns:
            C (numpy ndarray): The transformation matrix.
        """

        i =1
        j = 1
        k = 1
        if xi < (0-self.eps):
            i = -1
        if eta <(0-self.eps):
            j = -1
        if zeta <(0-self.eps):
            k = -1

        C = np.array([[i,0,0],[0,j,0],[0,0,k]])
            
        return C

    def _find_C4(self,xi,eta,zeta):
        """Finds the correct transformation matrix given the values of xi, eta, 
        and zeta for step 4.

        Args:
            xi (float): The value of xi.
            eta (float): The value of eta.
            zeta (float): The value of zeta.

        Returns:
            C (numpy ndarray): The transformation matrix.
        """

        i = 1
        j = 1
        k = 1
        p = None
        if (xi-self.eps)>0:
            i = -1
        elif not xi<(0-self.eps):
            p = 0
            
        if (eta-self.eps)>0:
            j = -1
        elif not eta<(0-self.eps):
            p = 1

        if (zeta-self.eps)>0:
            k = -1
        elif not zeta<(0-self.eps):
            p = 2

        if i*j*k<0:
            if p==0:
                i = -1
            elif p==1:
                j = -1
            else:
                k = -1

        C = np.array([[i,0,0],[0,j,0],[0,0,k]])

        return C   
