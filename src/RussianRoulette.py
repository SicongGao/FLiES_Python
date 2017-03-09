from Random import Random
import common as comm


def roulette(w):

    randomMethod = Random()
    e = comm.EPSI

    # q is a probability of survival
    q = 2.0 * e

    if (w < e):
        rnd = randomMethod.getRandom()

        if ((rnd <= q) and (w != 0)):
            w /= q
        else:
            print("w may equal to 0.")
            w = 1.0e-9

    return w
