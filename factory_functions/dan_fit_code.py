def fit_to_bare_model(data, yvar='t1'):
    def _adaptive_bare_model(x,p):
        k4 = x['k4']
        k6 = x['k6']
        beta = x['beta']

        y = {'a' : 0.}
        
        for param, value in p.items():
            if param[0] not in ['Q', 'C', 'I', 'L']:
                continue
            term = value
            if param[0] in ['Q', 'L']:
                for pp in param[1:]:
                    if pp == '4':
                        term *= k4
                    elif pp == '6':
                        term *= k6
                    elif pp == 'b':
                        term *= beta
            elif param[0] in ['I']:
                for pp in param[1:]:
                    if pp == '4':
                        term /= k4
                    elif pp == '6':
                        term /= k6
                    elif pp == 'b':
                        term /= beta
            y['a'] += term
        
        return y
    
    X = {
        'beta' : data.beta.values,
        'k4' : data.k4.values,
        'k6' : data.k6.values
    }
    # Switch to bare quark masses
    X['k4'] = 1./(2. * (4 - X['k4']))
    X['k6'] = 1./(2. * (4 - X['k6']))
    X['beta'] = 8. / X['beta']
    normalize_X(X)
    
    
    
    Y = {
        #'a' : 1./np.sqrt(data[yvar].values)
        #'a' : 1./data[yvar].values
        #'a' : data[yvar].values
        'a' : data[yvar].values
    }
    
    
    P = {
        'C' : gv.gvar(1,1e5),

        'Lb' : gv.gvar(1,1e5),
        'L4' : gv.gvar(1,1e5),
        'L6' : gv.gvar(1,1e5),

        'Qbb' : gv.gvar(1,1e5), 
        'Q44' : gv.gvar(1,1e5),
        'Q66' : gv.gvar(1,1e5),
        'Q4b' : gv.gvar(1,1e5),
        'Q6b' : gv.gvar(1,1e5),
        'Q46' : gv.gvar(1,1e5),
        
        'Qbbb' : gv.gvar(1,1e5),
        'Q444' : gv.gvar(1,1e5),
        'Q666' : gv.gvar(1,1e5),
        'Qbb4' : gv.gvar(1,1e5),
        'Qbb6' : gv.gvar(1,1e5),
        'Q446' : gv.gvar(1,1e5),
        'Q44b' : gv.gvar(1,1e5),
        'Q664' : gv.gvar(1,1e5),
        'Q66b' : gv.gvar(1,1e5),
        'Qb46' : gv.gvar(1,1e5),

        
        'Ib' : gv.gvar(1,1e5),
        'I4' : gv.gvar(1,1e5),
        'I6' : gv.gvar(1,1e5),

        'Ibb' : gv.gvar(1,1e5),
        'I44' : gv.gvar(1,1e5),
        'I66' : gv.gvar(1,1e5),

        'Ib4' : gv.gvar(1,1e5),
        'Ib6' : gv.gvar(1,1e5),
        'I46' : gv.gvar(1,1e5),
        
#         'Ibbb' : gv.gvar(1,1e5),
#         'I444' : gv.gvar(1,1e5),
#         'I666' : gv.gvar(1,1e5),
#         'Ibb4' : gv.gvar(1,1e5),
#         'Ibb6' : gv.gvar(1,1e5),
#         'I446' : gv.gvar(1,1e5),
#         'I44b' : gv.gvar(1,1e5),
#         'I664' : gv.gvar(1,1e5),
#         'I66b' : gv.gvar(1,1e5),
#         'Ib46' : gv.gvar(1,1e5),
    }
    
    
    cf = lsqfit.nonlinear_fit(data=(X,Y), fcn=_adaptive_bare_model,
#                               prior=P,
                              p0=gv.mean(P),
                              extend=True, debug=True)
    print cf.format(maxline=True)
    return cf