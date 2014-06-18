

def configuration(parent_package='', top_path=None):
    global config
    from numpy.distutils.misc_util import Configuration
    config = Configuration('retrieve', parent_package, top_path)
    config.add_data_dir('tests')

    print "FOOOOVAR"
    # XXX make this optional
    config.add_extension(
        'echo_steiner',
        sources=[generate_source],
        #sources=['echo_steiner.pyf', 'src/echo_steiner.f90'],
    )
    #print config.have_f90c()
    return config

def generate_source(ext, build_dir):
    print "KAHHN"
    if config.have_f90c():
        print "HAVE"
        return ['echo_steiner.pyf', 'src/echo_steiner.f90']
    else:
        print "HAVE NOT"
        return None

if __name__ == '__main__':
    from numpy.distutils.core import setup
    setup(**configuration(top_path='').todict())
