

def configuration(parent_package='', top_path=None):
    from numpy.distutils.misc_util import Configuration
    config = Configuration('retrieve', parent_package, top_path)
    #config.add_data_dir('tests')

    # echo classifier module
    config.add_extension(
        '_echo_classification',
        sources=['src/steiner.f90'])
    return config


if __name__ == '__main__':
    from numpy.distutils.core import setup
    setup(**configuration(top_path='').todict())
