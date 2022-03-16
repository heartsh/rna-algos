use utils::*;
pub const MATCH_SCORE_MAT: MatchScoreMat = [[0.5256508867, -0.40906402, -0.2502759109, -0.3252306723], [-0.40906402, 0.6665219366, -0.3289391181, -0.1326088918], [-0.2502759109, -0.3289391181, 0.6684676551, -0.3565888168], [-0.3252306723, -0.1326088918, -0.3565888168, 0.459052045]];
pub const INSERT_SCORES: InsertScores = [-0.002521927159, -0.08313891561, -0.07443970653, -0.01290054598];
pub const INIT_MATCH_SCORE: Prob = 0.3959924457;
pub const INIT_INSERT_SCORE: Prob = -0.3488104904;
pub const MATCH_2_MATCH_SCORE: Prob = 2.50575671;
pub const MATCH_2_INSERT_SCORE: Prob = 0.1970448791;
pub const INSERT_EXTEND_SCORE: Prob = 1.014026583;
pub const INSERT_SWITCH_SCORE: Prob = -7.346968782;