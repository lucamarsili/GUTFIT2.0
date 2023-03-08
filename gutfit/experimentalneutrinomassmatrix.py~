import numpy as np
from gutfit import model, parameterlist


#change predictions according to typeseesaw_v4 model

def matrix_diag3(d1,d2,d3):
    return np.array([[d1, 0.0, 0.0], [0.0, d2, 0.0], [0.0, 0.0, d3]])

# Generic Rotations #
def matrix_rot23(th23):
    s23, c23 = np.sin(th23), np.cos(th23)
    return np.array([[1,    0,   0],
                    [ 0,  c23, s23],
                    [ 0, -s23, c23]])

def matrix_rot12(th12):
    s12, c12 = np.sin(th12), np.cos(th12)
    return np.array([[ c12, s12, 0],
                    [ -s12, c12, 0],
                    [    0,   0, 1]])

def matrix_rot13(th13, delta):
    return np.array([[                     np.cos(th13), 0.0, np.sin(th13) * np.exp(-1j * delta)],
                    [                     0.0         , 1.0, 0.0                               ],
                    [-np.sin(th13)* np.exp(1j * delta), 0.0, np.cos(th13)]],
                    dtype=np.complex64)

def matrix_vckm(th12, th13, th23, delta):
    return matrix_rot23(th23) @ matrix_rot13(th13, delta) @ matrix_rot12(th12)


class ExperimentalNeutrinoMassMatrix(model.Model):
    def __init__(self):
        params = ["data_lepton_ye",
                  "data_lepton_ymu",
                  "data_lepton_ytau",
                  "data_neutrino_deltamsq21bf",
                  "data_neutrino_deltamsq31bf",
                  "data_neutrino_th23",
                  "data_neutrino_th12",
                  "data_neutrino_th13",
                  "data_neutrino_delta"
                  
                ]
        super().__init__(params)

    @property
    def val(self):
        value =  np.abs(
            self.Data(
                self.data_lepton_ye,
                self.data_lepton_ymu,
                self.data_lepton_ytau,
                self.data_neutrino_deltamsq21bf,
                self.data_neutrino_deltamsq31bf,
                self.data_neutrino_th23,
                self.data_neutrino_th12,
                self.data_neutrino_th13,
                self.data_neutrino_delta
            )
        )

        return value
    
    def MnuData(self, ye, ymu, ytau, deltamsq21bf,deltamsq31bf, th23l, th12l, th13l, deltal):
        m1logged = 10**m1
        mnudiag  = matrix_diag3(m1logged, np.sqrt(m1logged * m1logged+ deltamsq21bf), np.sqrt(m1logged * m1logged + deltamsq31bf))
        Majorana = matrix_diag3(1, 1,1)
        angle    = matrix_vckm(th12l, th13l, th23l, deltal)
        Vpmns    = angle 
        return np.conj(Vpmns) @ mnudiag @ np.transpose(Vpmns)
    
    
    def Data(self, ye, ymu, ytau,deltamsq21bf,deltamsq31bf, th23l, th12l, th13l, deltal):
        D = []
        D.append(th13l)
        D.append(th12l)
        D.append(th23l)
        D.append(deltamsq21bf)
        D.append(deltamsq31bf)
        D.append(ye)
        D.append(ymu)
        D.append(ytau)
        D = np.asarray(D)
        return D
       
       


if __name__=="__main__":
    E = ExperimentalNeutrinoMassMatrix()
    PL =  parameterlist.ParameterList.fromConfigFile("/home/lucamarsili/Documenti/Durham/gutfit-main/examples/param_card3sigma_v4.dat")
    #E(PL())
    import time
    t0 = time.time()
    print(E(PL()))
    print(time.time() - t0)
