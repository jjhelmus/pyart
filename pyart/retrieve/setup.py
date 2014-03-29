from os.path import join

def configuration(parent_package='', top_path=None):
    global config
    from numpy.distutils.misc_util import Configuration
    config = Configuration('retrieve', parent_package, top_path)
    config.add_data_dir('tests')

    # Conditionally add Steiner echo classifier extension.
    config.add_extension('_echo_steiner', sources=[steiner_echo_gen_source])
    config.add_extension('background',
                         sources=['background.pyf','src/background.f90'],
                         f2py_options=[])
    config.add_extension('continuity',
                         sources=['continuity.pyf','src/continuity.f90'],
                         f2py_options=[])
    config.add_extension('divergence',
                         sources=['divergence.pyf','src/divergence.f90'],
                         f2py_options=[])
    config.add_extension('gradient',
                         sources=['gradient.pyf','src/gradient.f90'],
                         f2py_options=[])
    config.add_extension('laplace',
                         sources=['laplace.pyf','src/laplace.f90'],
                         f2py_options=[])
    config.add_extension('smooth',
                         sources=['smooth.pyf','src/smooth.f90'],
                         f2py_options=[])

    return config


def steiner_echo_gen_source(ext, build_dir):
    """
    Add Steiner echo classifier source if Fortran 90 compliler available,
    if not compiler is found do not try to build the extension.
    """
    try:
        config.have_f90c()
        return [join(config.local_path, '_echo_steiner.pyf'),
                join(config.local_path, 'src', 'echo_steiner.f90')]
    except:
        # TODO add printer message about missing extension and
        # instruction on how to make available
        return None

if __name__ == '__main__':
    from numpy.distutils.core import setup
    setup(**configuration(top_path='').todict())
