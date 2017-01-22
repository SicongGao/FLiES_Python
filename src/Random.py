# ******************************************************************
#     Random number function based on the Tausworthe sequence method
#     Original code and algorithm:
#     Shinomoto (1992), Statistical mechanics: Parity physics cource.,
#     Maruzen press, P27
# *******************************************************************

# !!!!!!! HAVEN'T FINISH


class Random:

    r = 0.0
    A_NORM = 1.0 / 2147483647
    N_RAND = 1000
    cnt = 0

    # to make the sequence many times, following is repeated.
    def getRandom(self):

        if (self.cnt > self.N_RAND):
            for i in range(self.N_RAND):
                k = i * 100
                for j in range(i + k, k + 100):
                    return 0


        # remember to return a float
        return 0
