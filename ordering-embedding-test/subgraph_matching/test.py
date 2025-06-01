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

    # num_process = 20
    # process = []
    # Q = mp.Queue()
    # for i in range(num_process):
    #     p = mp.Process(target=load_query_graph, args=(model, test_pts[i][1], Q))
    #     p.start()
    #     process.append(p)

    # for p in process:
    #     p.join()

    # for i in range(num_process):
    #     emb_query = Q.get()
    #     if i == 0:
    #         emb_bs = emb_query
    #     else:
    #         emb_bs = torch.cat((emb_bs, emb_query), dim=0)

    
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
    

    # start_time = time.perf_counter()

    pred = model(emb_as, emb_bs)
    raw_pred = model.predict(pred)
    pred = model.clf_model(raw_pred.unsqueeze(1)).argmax(dim=-1)
    # raw_pred *= -1

    # print(raw_pred)

    end_time = time.perf_counter()
    confidence = torch.sigmoid(raw_pred)
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

    with open("./dataset_wordnet18/confidence.txt", "a") as fff:
        for i in range(78):
            fff.write(str(confidence[i].item()) + ' ')
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

    # pred = torch.cat(all_preds, dim=-1)
    # labels = torch.cat(all_labels, dim=-1)
    # raw_pred = torch.cat(all_raw_preds, dim=-1)
    # acc = torch.mean((pred == labels).type(torch.float))
    # prec = (torch.sum(pred * labels).item() / torch.sum(pred).item() if
    #     torch.sum(pred) > 0 else float("NaN"))
    # recall = (torch.sum(pred * labels).item() /
    #     torch.sum(labels).item() if torch.sum(labels) > 0 else
    #     float("NaN"))
    # labels = labels.detach().cpu().numpy()
    # raw_pred = raw_pred.detach().cpu().numpy()
    # pred = pred.detach().cpu().numpy()
    # auroc = roc_auc_score(labels, raw_pred)
    # avg_prec = average_precision_score(labels, raw_pred)
    # tn, fp, fn, tp = confusion_matrix(labels, pred).ravel()
    # if verbose:
    #     import matplotlib
    #     matplotlib.use('Agg')
    #     import matplotlib.pyplot as plt
    #     precs, recalls, threshs = precision_recall_curve(labels, raw_pred)
    #     plt.plot(recalls, precs)
    #     plt.xlabel("Recall")
    #     plt.ylabel("Precision")
    #     plt.savefig("plots/precision-recall-curve.png")
    #     print("Saved PR curve plot in plots/precision-recall-curve.png")

    print("\n{}".format(str(datetime.now())))
    # print("Validation. Epoch {}. acc: {:.4f}. "
    #     "pre: {:.4f}. recall: {:.4f}. f1: {:.4f}.\n     "
    #     "TN: {}. FP: {}. FN: {}. TP: {}".format(epoch,
    #         acc, pre, recall, f1,
    #         tn, fp, fn, tp))

    # if not args.test:
    #     logger.add_scalar("Accuracy/test", acc, batch_n)
    #     logger.add_scalar("Precision/test", prec, batch_n)
    #     logger.add_scalar("Recall/test", recall, batch_n)
    #     logger.add_scalar("AUROC/test", auroc, batch_n)
    #     logger.add_scalar("AvgPrec/test", avg_prec, batch_n)
    #     logger.add_scalar("TP/test", tp, batch_n)
    #     logger.add_scalar("TN/test", tn, batch_n)
    #     logger.add_scalar("FP/test", fp, batch_n)
    #     logger.add_scalar("FN/test", fn, batch_n)
    #     print("Saving {}".format(args.model_path))
    #     torch.save(model.state_dict(), args.model_path)

    # if verbose:
    #     conf_mat_examples = defaultdict(list)
    #     idx = 0
    #     for pos_a, pos_b, neg_a, neg_b in test_pts:
    #         if pos_a:
    #             pos_a = pos_a.to(utils.get_device())
    #             pos_b = pos_b.to(utils.get_device())
    #         neg_a = neg_a.to(utils.get_device())
    #         neg_b = neg_b.to(utils.get_device())
    #         for list_a, list_b in [(pos_a, pos_b), (neg_a, neg_b)]:
    #             if not list_a: continue
    #             for a, b in zip(list_a.G, list_b.G):
    #                 correct = pred[idx] == labels[idx]
    #                 conf_mat_examples[correct, pred[idx]].append((a, b))
    #                 idx += 1

if __name__ == "__main__":
    from subgraph_matching.train import main
    main(force_test=True)
