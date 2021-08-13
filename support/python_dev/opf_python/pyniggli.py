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

    def __init__(self,A,eps = None, path_=None):
        """Initial setup of cell.
        Args:
            A (numpy ndarray): A 3 by 3 matrix containing the lattice vectors as columns.
            eps (optional float): Floating point tollerance for comparisons, 
                default is 1E-5.
            path_ (optional bool): True if the path of the reduction should be printed
                after the reduction, default False.

        Rasise:
            ValueError: if the input is not a 3 by 3 matrix.
            ValueError: if the input has a determinant of 0.
            RuntimeError: if the niggli cell is not found within 100 iterations.
        """

        if path_ is None:
            path = False
        else:
            path = path_
            
        if not isinstance(A,np.ndarray):
            A = np.array(A)
        if A.shape != (3,3):
            raise ValueError("The input basis must be a 3 by 3 matrix with "
                             "the lattice vectors as columns of the matrix.")
        
        if np.linalg.det(A) ==0:
            raise ValueError("The cell specified has a volume of zero.")
        else:
            self.volume = abs(np.linalg.det(A))
            
        self.original = np.array(A)
        if eps is None:
            self.eps = (1E-5)*self.volume**(1./3.)
        else:
            self.eps = eps*self.volume**(1./3.)
        self._niggli_reduction(path)
        self.niggli = np.dot(self.original,self.C)


    def _get_params(self):
        """Gets the niggli parameters A, B, C, xi, eta, zeta, l, m, n.

        Returns:
            A, B, C, xi, eta, zeta, l, m, n (float x 6, int x 3): The niggli
                parameters.
        """

        mat = np.dot(self.original,self.C)
        print("mat",mat)
        a = mat[:,0]
        b = mat[:,1]
        c = mat[:,2]

        A = np.dot(a,a)
        B = np.dot(b,b)
        C = np.dot(c,c)
        xi = 2.0*np.dot(b,c)
        eta = 2.0*np.dot(c,a)
        zeta = 2.0*np.dot(a,b)
        l, m, n = 0, 0, 0
        if xi<-self.eps:
            l = -1
        elif xi> self.eps:
            l = 1
            
        if eta<-self.eps:
            m = -1
        elif eta> self.eps:
            m = 1
                
        if zeta<-self.eps:
            n = -1
        elif zeta> self.eps:
            n = 1

        return A, B, C, xi, eta, zeta, l, m, n
        
    def _niggli_reduction(self,print_path):
        """Performs the niggli reduction of the given lattice.

        Args:
            path (bool): True if the path of the reduction should be printed.

        Raises:
            RuntimeError: if the niggli cell is not found within 100 iterations.
        """

        count, reduced, self.C = 0, False, np.array([[1,0,0],[0,1,0],[0,0,1]]) 

        A, B, C, xi, eta, zeta, l, m, n = self._get_params()

        path = ""
        while not reduced and count <=1000:
                
            reduced = True
            count += 1
            #1
            if A > (B+self.eps) or (not (abs(A-B)>self.eps) and abs(xi)>(abs(eta)+self.eps)):
                path += "1"
                self.C = np.dot(self.C,[[0,-1,0],[-1,0,0],[0,0,-1]])
                A, B, C, xi, eta, zeta, l, m, n = self._get_params()
                reduced = False
                
            #2
            if B > (C+self.eps) or (not (abs(C-B)>self.eps) and abs(eta)>abs(zeta)+self.eps):
                path += "2"
                self.C = np.dot(self.C,[[-1,0,0],[0,0,-1],[0,-1,0]])
                A, B, C, xi, eta, zeta, l, m, n = self._get_params()
                reduced = False
                continue
                #go to 1
            #3
            if l*m*n==1:
                path += "3"
                M = self._find_C3(l,m,n)
                self.C = np.dot(self.C,M)
                A, B, C, xi, eta, zeta, l, m, n = self._get_params()
                if not np.allclose(M,[[1,0,0],[0,1,0],[0,0,1]]):
                    reduced = False
                
            #4
            if l*m*n == 0 or l*m*n == (-1):
                path += "4"
                M = self._find_C4(l,m,n)
                self.C = np.dot(self.C,M)
                if not np.allclose(M,[[1,0,0],[0,1,0],[0,0,1]]):
                    reduced = False
                A, B, C, xi, eta, zeta, l, m, n = self._get_params()
                
            #5
            if abs(xi)>(B+self.eps) or (not (abs(B-xi)>self.eps) and (2*eta<(zeta-self.eps))) or (not abs(B+xi)>self.eps and zeta<(-self.eps)):
                path += "5"
                self.C = np.dot(self.C,np.array([[1,0,0],[0,1,-np.sign(xi)],[0,0,1]]))
                A, B, C, xi, eta, zeta, l, m, n = self._get_params()
                reduced = False
                continue
                #go to 1
            #6
            if abs(eta)>(A+self.eps) or (not abs(A-eta)>self.eps and (2*xi<(zeta-self.eps))) or (not abs(A+eta)>self.eps and zeta<(-self.eps)):
                path += "6"
                self.C = np.dot(self.C,np.array([[1,0,-np.sign(eta)],[0,1,0],[0,0,1]]))
                A, B, C, xi, eta, zeta, l, m, n = self._get_params()
                reduced = False
                continue
                #go to 1
            #7
            if abs(zeta)>(A+self.eps) or (not abs(A-zeta)>self.eps and (2*xi<(eta-self.eps))) or (not abs(A+zeta)>self.eps and eta<(-self.eps)):
                path += "7"
                self.C = np.dot(self.C,np.array([[1,-np.sign(zeta),0],[0,1,0],[0,0,1]]))
                A, B, C, xi, eta, zeta, l, m, n = self._get_params()
                reduced = False
                continue
                #go to 1
            #8
            if xi+eta+zeta+A+B<(-self.eps) or (not abs(xi+eta+zeta+A+B)>self.eps and (2*(A+eta)+zeta)>self.eps):
                path += "8"
                self.C = np.dot(self.C,np.array([[1,0,1],[0,1,1],[0,0,1]]))
                A, B, C, xi, eta, zeta, l, m, n = self._get_params()
                reduced = False
                continue
                #go to 1

        if print_path:
            print(path)

        if count >= 1000: #pragma: no cover
            raise RuntimeError("Could not reduce the cell in 1000 iterations. This could be "
                               "because of floating point error, try providing a smaller eps "
                               "value.", self.original)

        if not self._niggli_check(A,B,C,xi,eta,zeta,self.eps): # pragma: no cover
            raise RuntimeError("Cell reduction incorroct A: {0}, B: {1}, C: {2}, "
                               "xi: {3}, eta: {4}, zeta: {5}. Please submit a bug "
                               "report to: https://github.com/wsmorgan/pyniggli/"
                               "issues".format(A,B,C,xi,eta,zeta), self.original)


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
            return False
        
        if np.allclose(A,B,atol=eps) and not (abs(xi) < abs(eta)-eps or
                                              np.allclose(abs(xi),abs(eta),atol=eps)):
            return False
            
        if np.allclose(B,C,atol=eps) and not (abs(eta) < abs(zeta)-eps
                                              or np.allclose(abs(eta),abs(zeta),atol=eps)):
            return False

        if not ((xi-eps > 0 and eta-eps > 0 and zeta-eps > 0) or
                ((xi < 0-eps or np.allclose(xi,0,atol=eps))
                 and (eta < 0-eps or np.allclose(eta,0,atol=eps))
                 and (zeta < 0-eps or np.allclose(zeta,0,atol=eps)))):
            return False

        if not (abs(xi) < B-eps or np.allclose(abs(xi),B,atol=eps)):
            return False 
            
        if not ((abs(eta) < A-eps or np.allclose(abs(eta),A,atol=eps)) and (abs(zeta) < A-eps or
                                                               np.allclose(abs(zeta),A,atol=eps))):
            return False

        if not (C < A+B+C+xi+eta+zeta-eps or np.allclose(C, A+B+C+xi+eta+zeta,atol=eps)):
            return False

        if np.allclose(xi,B,atol=eps) and not (zeta < 2.*eta-eps or
                                               np.allclose(zeta,2.*eta,atol=eps)):
            return False

        if np.allclose(eta,A,atol=eps) and not (zeta < 2.*xi-eps or
                                                np.allclose(zeta,2.*xi,atol=eps)):
            return False

        if np.allclose(zeta,A,atol=eps) and not (eta < 2.*xi-eps or
                                                 np.allclose(eta,2.*xi,atol=eps)):
            return False

        if np.allclose(xi,-B,atol=eps) and not np.allclose(zeta,0,atol=eps):
            return False

        if np.allclose(eta,-A,atol=eps) and not np.allclose(zeta,0,atol=eps):
            return False

        if np.allclose(zeta,-A,atol=eps) and not np.allclose(eta,0,atol=eps):
            return False

        if np.allclose(C,A+B+C+xi+eta+zeta,rtol=0.0) and not ((2.*A+2.*eta+zeta) < 0-eps or
                                                     np.allclose(2.*A+2.*eta+zeta,0,atol=eps)):
            return False

        return True
    
    @staticmethod
    def _find_C3(l,m,n):
        """Finds the correct transformation matrix given the values of xi, eta, 
        and zeta for step 3.

        Args:
            l (int): Positive if xi is positive, negative if it isn't.
            m (int): Positive if eta is positive, negative if it isn't.
            n (int): Positive if zeta is positive, negative if it isn't.

        Returns:
            C (numpy ndarray): The transformation matrix.
        """

        i =1
        j = 1
        k = 1
        if l == (-1):
            i = -1
        if m == (-1):
            j = -1
        if n == (-1):
            k = -1

        C = np.array([[i,0,0],[0,j,0],[0,0,k]])
            
        return C

    @staticmethod
    def _find_C4(l,m,n):
        """Finds the correct transformation matrix given the values of xi, eta, 
        and zeta for step 4.

        Args:
            l (int): Positive if xi is positive, negative if it isn't.
            m (int): Positive if eta is positive, negative if it isn't.
            n (int): Positive if zeta is positive, negative if it isn't.

        Returns:
            C (numpy ndarray): The transformation matrix.
        """

        i = 1
        j = 1
        k = 1

        if l == (-1) and m == (-1) and n == (-1):
            C = np.array([[i,0,0],[0,j,0],[0,0,k]])
            return C
        
        r = -1
        if l == 1:
            i = -1
        elif l == 0:
            r = 0

        if m == 1:
            j = -1
        elif m == 0:
            r = 1

        if n == 1:
            k = -1
        elif n == 0:
            r = 2

        if i*j*k == (-1):
            if r==0:
                i = -1
            elif r==1:
                j = -1
            elif r==2:
                k = -1
                
        C = np.array([[i,0,0],[0,j,0],[0,0,k]])

        return C   
