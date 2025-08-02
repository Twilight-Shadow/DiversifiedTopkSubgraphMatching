from common import utils
from collections import defaultdict
from datetime import datetime
from sklearn.metrics import roc_auc_score, confusion_matrix
from sklearn.metrics import precision_recall_curve, average_precision_score
import torch
import torch.multiprocessing as mp
import time
#matplotlib.use('Agg')

USE_ORCA_FEATS = False # whether to use orca motif counts along with embeddings
MAX_MARGIN_SCORE = 1e9 # a very large margin score to given orca constraints

def load_query_graph(model, test_pts, queue):
    test_pts = test_pts.to(utils.get_device())
    with torch.no_grad():
        emb_query = model.emb_model(test_pts)
    queue.put(emb_query)

def validation(args, model, test_pts, logger, batch_n, epoch, verbose=False):
    # test on new motifs
    model.eval()
    all_raw_preds, all_preds = [], []

    print(len(test_pts))

    for i in range(len(test_pts)):
        neg_a = test_pts[i][2]
        neg_a = neg_a.to(utils.get_device())
        with torch.no_grad():
            emb_data = model.emb_model(neg_a)
            if i == 0:
                emb_as = emb_data
            else:
                emb_as = torch.cat((emb_as, emb_data), dim=0)

    start_time = time.perf_counter()
    
    for i in range(1):
        neg_b = test_pts[i][1]
        neg_b = neg_b.to(utils.get_device())
        with torch.no_grad():
            emb_query = model.emb_model(neg_b)
            if i == 0:
                emb_bs = emb_query
            else:
                emb_bs = torch.cat((emb_bs, emb_query), dim=0)

    for i in range(0):
        emb_bs = torch.cat((emb_bs, emb_query), dim=0)

    pred = model(emb_as, emb_bs)
    raw_pred = model.predict(pred)
    pred = model.clf_model(raw_pred.unsqueeze(1)).argmax(dim=-1)

    end_time = time.perf_counter()
    elapsed_time = end_time - start_time
    print(f"time cost: {elapsed_time} s")

    tp = tn = fp = fn = 0
    ground_truth = test_pts[0][4]
    # print(len(ground_truth))
    # for i in range(78, 1600):
    #     if ground_truth[i] == 1 and pred[i - 78] == 1:
    #         tp += 1
    #     if ground_truth[i] == 0 and pred[i - 78] == 1:
    #         fp += 1
    #     if ground_truth[i] == 1 and pred[i - 78] == 0:
    #         fn += 1
    #     if ground_truth[i] == 0 and pred[i - 78] == 0:
    #         tn += 1

    # print(ground_truth[:1000])
    # print(pred)
    # print(pred.shape)

    with open("./dataset_wordnet18/pred.txt", "a") as fff:
        for i in range(78):
            fff.write(str(pred[i].item()) + ' ')
        fff.write("\n")

    if tp + fp != 0:
        pre = tp / (tp + fp)
    else:
        pre = float("NaN")
    if tp + fn != 0:
        recall = tp / (tp + fn)
    else:
        recall = float("NaN")
    acc = (tp + tn) / 78
    if pre != float("NaN") and recall != float("NaN"):
        f1 = (2 * pre * recall) / (pre + recall)
    else:
        f1 = float("NaN")

    print("\n{}".format(str(datetime.now())))

if __name__ == "__main__":
    from subgraph_matching.train import main
    main(force_test=True)
