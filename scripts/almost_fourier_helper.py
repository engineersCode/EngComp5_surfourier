import numpy
from IPython.display import HTML
import ipywidgets
from matplotlib import animation, pyplot

def create_init_fig(wrapped_signal, freq_arr, xcm_arr):
    """ creates initial figure needed for animation, but it doesn't display it.
    """
    
    pyplot.ioff() #to avoid displaying figure. 

    fig, ax = pyplot.subplots(figsize=(14.0, 6.0))
    pyplot.tight_layout()
    fig.suptitle('Frequency = {:.2f}'.format(freq_arr[0]))

    ax1 = pyplot.subplot2grid((1, 3), (0, 0))
    ax2 = pyplot.subplot2grid((1, 3), (0, 1), colspan=2)

    circle1 = pyplot.Circle((0, 0), 1, fill=None, lw=2, ls='--', alpha=0.3)

    ax1.add_patch(circle1)
    ax1.grid()

    ticks= numpy.linspace(-1,1, 5, endpoint=True)

    ylabels = [-1, -0.5, None, 0.5, 1]

    ax1.set_yticklabels(ylabels)
    ax1.set_xticks(ticks)
    ax1.set_yticks(ticks)

    wrapped_signal_plot = ax1.plot(wrapped_signal.real, 
                                   wrapped_signal.imag, alpha=0.5,
                                   label=r'$g(t)e^{2\pi ift}$')[0]

    # Move left y-axis and bottim x-axis to centre, passing through (0,0)
    ax1.spines['left'].set_position('center')
    ax1.spines['bottom'].set_position('center')

    # Eliminate upper and right axes
    ax1.spines['right'].set_color('none')
    ax1.spines['top'].set_color('none')


    ax1.axis('scaled', adjustable='box')
    ax1.set_xlim(-1.1,1.1)
    ax1.set_ylim(-1.1,1.1)
    ax1.legend(loc='upper left', bbox_to_anchor=(0.48, 1.12))

    #f_list = numpy.full_like(freqs, None)
    almost_fourier_plot = ax2.plot(freq_arr[0],  xcm_arr[0], '-')[0]
    ax2.spines['right'].set_color('none')
    ax2.spines['top'].set_color('none')
    ax2.axis('scaled', adjustable='box')
    ax2.set_xlabel('Frequency')
    ax2.set_ylabel('xcm', labelpad=-15)

    ax2.set_xlim(0.9,4.1)
    ax2.set_ylim(-0.3,1.1)
    ax2.grid();
    
    return {'fig': fig, 'WSP': wrapped_signal_plot, 'AF': almost_fourier_plot}


def update_figure(f, anim_dict, g_t, freq_arr, t_arr, display_fig=False):
    
    res = g_t * anim_dict['_comp_fourier'](freq_arr[f], t_arr)
    
    anim_dict['fig'].suptitle('Frequency = {:.2f}'.format(freq_arr[f]))
    
    anim_dict['xcm_arr'][f] = anim_dict['_get_xcm'](res)
    
    anim_dict['WSP'].set_data(res.real, res.imag)
    
    anim_dict['AF'].set_data(freq_arr[:f+1], anim_dict['xcm_arr'][:f+1])
    
    if display_fig:
        display(anim_dict['fig'])