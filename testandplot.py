import pickle
import sys




def main():
    data = pickle.load(file(sys.argv[1]))


if __name__ == "__main__":
    main()