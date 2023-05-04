import numpy as np
def rfbandwidth(Z,nrf,loss=1.0,lengthinmm=10,λ=1064,nktp=1.8,fmax=40,df=0.01,za=50,zt=50,legendtext='',plot=True): # loss in dB/cm/√GHz, frequnecy in GHz, λ in nm, Z in Ω, nktp = refractive index
    # 20*log for optical response, 20*log for S11,S21 and RF loss, 10*log for |S11|²
    length = lengthinmm/1000
    x = np.arange(0,fmax+df/2,df)
    x[0] += 1e-9
    def α(index):
        return np.log(10**(loss/20))*np.sqrt(x)*length/.01 + 1j*length*2*np.pi*1e9*x/299792458*index
    def rfprop(za,z0,zt,alpha0): # za = input impedance, z0 = electrode impedance, zt = termination impedance
        g = (zt-z0)/(zt+z0)
        ge = np.exp(-alpha0)
        ga = ( 1/ge-g*ge )/( 1/ge+g*ge ) * za/z0
        ga = (1-ga)/(1+ga)
        gb = (1+g)*(1+ga)/( 1/ge+g*ge )
        # ga[0] = nan
        return g,ge,ga,gb
    def opprop(g,ge,ga,alpha1,alpha2):
        gg = (np.exp(+alpha1+1e-99)-1)/(+alpha1+1e-99)
        gg += g*(np.exp(-alpha2)-1)/(-alpha2)
        gg *= (1+ga)/( 1/ge+g*ge )
        return gg
    def dbwave(y,name=''):
        # return Wave(20*np.log10(abs(y)),x,name)
        return 20*np.log10(abs(y))
    za, z0, zt = za+0j, Z+0j, zt+0j
    g,ge,ga,gb = rfprop(za,z0,zt,α(nrf))
    opticalcoprop = dbwave( opprop(g,ge,ga,α(nrf-nktp),α(nrf+nktp)),'optical, co-prop' )
    opticalcounterprop = dbwave( opprop(g,ge,ga,α(nrf+nktp),α(nrf-nktp)),'optical, counter-prop' )
    opticalvelocitymatched = dbwave( opprop(g,ge,ga,α(0),α(2*nrf)),'optical, velocity matched' )
    uu = dbwave( np.sqrt(ga*ga.conjugate() + gb*gb.conjugate()), 'RF, |S11|²+|S21|²' )
    s11 = dbwave(ga,'RF, |S11|')
    s21 = dbwave(gb,'RF, |S21|')
    opticallossless = dbwave( opprop( *rfprop(za,z0,zt,1j*α(nrf).imag)[:3], 1j*α(nrf-nktp).imag, 1j*α(nrf+nktp).imag),'optical, no loss' )
    def sinc(x):
        return np.sinc(x/np.pi)
    def bestlengthinmm(f):
        return 1000*299792458/(2*1e9*f*(nrf-nktp))
    opticallossless50ohm = dbwave( sinc(0.5*α(nrf-nktp).imag),'optical, no loss 50Ω' )
    opticallossless50ohmcounterprop = dbwave( sinc(0.5*α(nrf+nktp).imag),'optical, counter-prop no loss 50Ω' )
    opbestlen = dbwave( np.where(0.5*α(nrf-nktp).imag>np.pi/2, 1/(0.5*α(nrf-nktp).imag+1e-99), sinc(0.5*α(nrf-nktp).imag)),'optical, best length' ) # opticallossless50ohm with best length chosen, length = min( len, bestlength(f) )
    # if plot: Wave.plots(opticalcoprop,opticalcounterprop,opticalvelocitymatched,opticallossless,
    #         opticallossless50ohm,opticallossless50ohmcounterprop,opbestlen,s11,s21,uu,
    #         seed=16,x='frequency (GHz)',y='response (dB)',ylim=((np.nanmin(opticalcounterprop)//10)*10,None),
    #         legendtext=(legendtext+'\n' if legendtext else '')+f'RF loss = {loss}dB/cm/√GHz',
    #         save='rfbandwidth plot')
    # return {'coprop':opticalcoprop,'counterprop':opticalcounterprop,'velocitymatched':opticalvelocitymatched,'s11':s11,'s21':s21,'lossless':opticallossless,'lossless50ohm':opticallossless50ohm,'lossless50ohmcounterprop':opticallossless50ohmcounterprop,'opbestlen':opbestlen}
    return {'frequency':x,'coprop':opticalcoprop,'counterprop':opticalcounterprop,'velocitymatched':opticalvelocitymatched,'s11':s11,'s21':s21,'lossless':opticallossless,'lossless50ohm':opticallossless50ohm,'lossless50ohmcounterprop':opticallossless50ohmcounterprop,'opbestlen':opbestlen}
def plotit(d):
    import matplotlib.pyplot as plt
    # plt.plot(d['frequency'],d['coprop'])
    for s in d:
        if not s=='frequency':
            plt.plot(d['frequency'],d[s],label=s)
    plt.xlabel(u'frequency (GHz)')
    plt.ylabel(u'response (dB)')
    plt.legend()
    plt.show()
if __name__ == '__main__':
    d = rfbandwidth(Z=50,nrf=2.8,loss=1,za=40,zt=50,lengthinmm=10,λ=1064,nktp=1.8,fmax=40,plot=0)
    print(d['coprop'])
    plotit(d)
