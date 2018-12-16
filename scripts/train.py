import copy
import sklearn.metrics

def fitAndEval(X_train, Y_train, X_test, Y_test, model, copyModel=True, returnModel=False):
    '''
    Args
    - X_train, Y_train: array-like. shape=(n_samples, n_features)
    - Y_train, Y_test: list-like. shape=(n_samples,)
    - model: object
        Must have .fit() and decision_function() or predict_proba() methods
    - copyModel: bool. default=True
        False: Use `model` by reference; the fitted model is then accessible in the calling envrionment
        True: Make a deep copy of `model` before beginning. Fitting does not affect the model in the calling environment
    - returnModel: bool. default=False
        Include the fitted model as the 5th element of the returned tuple
    
    Returns: tuple
      returnModel==False: (fpr, tpr, thresh, roc_auc)
      returnModel==True: (fpr, tpr, thresh, roc_auc, model)
    '''
    
    if copyModel:
        model = copy.deepcopy(model)
    
    model.fit(X_train, Y_train)
    print("done fitting" + str(model))
    
    # get target scores to feed into sklearn.metrics.roc_curve()
    if hasattr(model, "decision_function") and callable(model.decision_function):
        # decision function (distance to separating hyperplane) for SVM models
        y_score_train = model.decision_function(X_train)
        y_score = model.decision_function(X_test)
    elif hasattr(model, "predict_proba") and callable(model.predict_proba):
        # probability estimates returned available for most sklearn models (linear, ensemble, MLPClassifier, etc.)
        y_score_train = model.predict_proba(X_train)[:,1]
        y_score = model.predict_proba(X_test)[:,1]
    
    fpr_train, tpr_train, thresh_train = sklearn.metrics.roc_curve(Y_train, y_score_train)
    roc_auc_train = sklearn.metrics.auc(fpr_train, tpr_train)
    fpr, tpr, thresh = sklearn.metrics.roc_curve(Y_test, y_score)
    roc_auc = sklearn.metrics.auc(fpr, tpr)
    if returnModel:
        return (fpr, tpr, thresh, roc_auc, fpr_train, tpr_train, thresh_train, roc_auc_train, model)
    else:
        return (fpr, tpr, thresh, roc_auc, fpr_train, tpr_train, thresh_train, roc_auc_train)

def fit(model, X, Y, copyModel=True, returnModel=True, verbose=True):
    if copyModel:
        model = copy.deepcopy(model)
    
    model.fit(X, Y)
    if verbose:
        print("done fitting" + str(model))
    
    if returnModel:
        return model

def eval(model, X, Y):
    '''
    Args
    - model
    - X
        dataset to predict on
    - Y
        labels of dataset to predict on
    '''
    
    # get target scores to feed into sklearn.metrics.roc_curve()
    if hasattr(model, "decision_function") and callable(model.decision_function):
        # decision function (distance to separating hyperplane) for SVM models
        y_score = model.decision_function(X)
    elif hasattr(model, "predict_proba") and callable(model.predict_proba):
        # probability estimates returned available for most sklearn models (linear, ensemble, MLPClassifier, etc.)
        y_score = model.predict_proba(X)[:,1]
    
    fpr, tpr, thresh = sklearn.metrics.roc_curve(Y, y_score)
    roc_auc = sklearn.metrics.auc(fpr, tpr)
    return {'fpr': fpr, 'tpr': tpr, 'thresh': thresh, 'roc_auc': roc_auc}