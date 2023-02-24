from argparse import ArgumentParser

download_argument_parser = ArgumentParser()
download_argument_parser.add_argument('directory')

preprocess_argument_parser = ArgumentParser()
preprocess_argument_parser.add_argument('directory')

train_model_argument_parser = ArgumentParser()
train_model_argument_parser.add_argument('directory')

search_argument_parser = ArgumentParser()
search_argument_parser.add_argument('-i', '--input', required=True)
search_argument_parser.add_argument('-m', '--model', required=True)
search_argument_parser.add_argument('-o', '--output', required=True)
search_argument_parser.add_argument('-f1', '--minimum-f1-score', default=0.65, type=float)
