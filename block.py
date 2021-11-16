class Block:
    def __init__(self):
        self.start = 0
        self.end = 0
        self.het_idxes = []
        self.het_poses = []
        self.whole_size = self.end-self.start
        self.het_size = len(self.het_idxes)
        self.het_seq = ""
        self.whole_seq = ""
        self.het_pos_seq = ""
        self.het_neg_seq = ""