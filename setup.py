from distutils.core import setup

setup(name = 'IGCexpansion',
      version = '0.3-dev',
      author = 'Xiang Ji',
      url = 'https://github.com/xji3/IGCexpansion',
      #download_url = 'https://github.com/xji3/Genconv/tree/master/IGCexpansion/',
      packages = ['IGCexpansion',],
      install_requires=[
          'Biopython', 'networkx', 'numpy', 'scipy', #'jsonctmctree'
      ],
      dependency_links=[
      'https://github.com/xji3/jsonctmctree/tarball/master#egg=package-0.2.0',
      'git+https://github.com/xji3/jsonctmctree.git@master#egg=jsonctmctree-0.2.0'
      ]
      #long_description = open('README.md').read()
      )
