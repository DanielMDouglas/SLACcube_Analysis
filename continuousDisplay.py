import matplotlib.pyplot as plt
import trackDisplay

if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument('infile',
                        type = str,
                        help = 'flattened hit data from selection.py')
    parser.add_argument('-t', '--trackID',
                        type = int,
                        default = 0,
                        help = 'trackID to plot')
    parser.add_argument('--plotfile', '-p',
                        default = "",
                        type = str,
                        help = 'optional file to which resulting track display is saved')
    args = parser.parse_args()

    while True:
        status = trackDisplay.main(args)
        if( status != 0 ):
            break
        key_input = input()
        if (key_input == 'k') or (key_input == 'n'):
            plt.close()
            args.trackID += 1
        elif (key_input == 'j') or (key_input == 'b'):
            plt.close()
            args.trackID -= 1
        elif (key_input == 'e') or (key_input == 'q'):
            plt.close()
            break
        else :
            print("Type \'n\' (\'b\') for next (previous) event, \'e\' or \'q\' to exit")
 