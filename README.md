# radar-tools
AMISR tools for radar related calculations. e.g. routines to calculate ambiguity function.

# install
Clone this repo and:

    $ pip install .
    
# Usage

## io_utils
Once installed in a python script or notebook:

    >>> from radartools import io_utils
    >>> f0 = io_utils.read_whole_h5file("acfl16_30_10_ph/AmbFunc.h5")
    >>> f0.print_info()
    
    Default argument: plotattrs=False
    /
      |__ Delay ............. float64...... (990,)..............     7.73 kB [0]=-0.000254 [-1]=0.000735
      |__ Lags .............. float64...... (48,)...............      384 Bytes [0]=0 [-1]=0.00047
      |__ Range ............. float64...... (735,)..............     5.74 kB [0]=0 [-1]=110024
      |__ Wlag .............. complex128... (48, 990)...........      742 kB
      |__ WlagSum ........... complex128... (48,)...............      768 Bytes [0]=(0.290413+-0.0217599j) [-1]=(0.0196897+-0.00377615j)
      |__ Wrange ............ complex128... (48, 735)...........      551 kB
      |__ WrangeSum ......... complex128... (48,)...............      768 Bytes [0]=(0.290413+-0.0217599j) [-1]=(0.0196897+-0.00377615j)
      
## amb_func_tools
See example in examples/acfl16_30_10_ph on how to generate AmbFunc.h5 for alternating codes with a phase correction file.
Additional plotting routines to plot any AmbFunc.h5 :

    >>> from radartools import io_utils
    >>> from radartools import amb_func_tools
    >>> f0 = io_utils.read_whole_h5file("acfl16_30_10_ph/AmbFunc.h5")
    >>> amb_func_tools.plot_complex_Wlag_Wrange(f0, title="acfl16_30_10_ph",
            figname="acfl16_30_10_ph_Wlag_Wrange.png")

