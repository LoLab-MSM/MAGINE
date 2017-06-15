

def translate(X, Y, option, param_counter, parameters, rules, gene_monomers,
              compound_monomers):
    """
    translate edge information to a rule


    Parameters
    ----------
    X : str
    Y : str
    option : str
    param_counter : int
    parameters : list
    rules : list
    gene_monomers : list
    compound_monomers : list

    Returns
    -------

    """
    X = str(X)
    Y = str(Y)
    split = option.split("_")

    if split[0] == 'inhibition':
        rate1 = "kf%s" % str(param_counter)
        rate2 = "kr%s" % str(param_counter)
        rate3 = "kc%s" % str(param_counter)
        if X in gene_monomers:
            enzyme = '%s(gene="protein",state="A")' % X
            sub = '%s(gene="protein",state="A")' % Y
            prod = '%s(gene="protein",state="I")' % Y
        else:
            enzyme = '%s(bf=None)' % X
            sub = '%s(gene="protein",state="A")' % Y
            prod = '%s(gene="protein",state="I")' % Y

        if Y in compound_monomers:
            enzyme = '%s(bf=None)' % Y
            sub = '%s(gene="protein",state="A")' % X
            prod = '%s(gene="protein",state="I")' % X

        tmp1 = 'catalyze(%s,%s,%s,[%s,%s,%s] )\n' % (
            enzyme, sub, prod, rate1, rate2, rate3)
        rules += tmp1
        parameters += 'Parameter("%s",1)\n' % rate1
        parameters += 'Parameter("%s",1)\n' % rate2
        parameters += 'Parameter("%s",1)\n' % rate3

    elif split[0] == 'expression':

        rate1 = "kf%s" % str(param_counter)
        rate2 = "kr%s" % str(param_counter)
        rate3 = "kc%s" % str(param_counter)
        rate4 = "k_synth%s" % str(param_counter)
        rule_name = '"%s_expression_%s"' % (Y, str(param_counter))
        tmp1 = 'catalyze(%s(gene="protein"),%s(gene="off"),%s(gene="on",state="A"),[%s,%s,%s] )\n' % (
            X, Y, Y, rate1, rate2, rate3)
        tmp2 = 'Rule(%s,%s(bf=None,gene="on",state="A") !> %s(bf=None,gene="on",state="A") + %s(bf=None,gene="protein",state="A"),%s)\n' % (
            rule_name, Y, Y, Y, rate2)
        rules += tmp1
        rules += tmp2
        parameters += 'Parameter("%s", 1)\n' % rate1
        parameters += 'Parameter("%s", 1)\n' % rate2
        parameters += 'Parameter("%s", 1)\n' % rate3
        parameters += 'Parameter("%s", 1)\n' % rate4

    elif split[0] == 'activation':

        rate1 = "kf%s" % str(param_counter)
        rate2 = "kr%s" % str(param_counter)
        rate3 = "kc%s" % str(param_counter)

        if X in gene_monomers:
            enzyme = '%s(gene="protein", state="A")' % X
            sub = '%s(gene="protein", state="I")' % Y
            prod = '%s(gene="protein", state="A")' % Y
        else:

            enzyme = '%s(bf=None)' % X
            sub = '%s(gene="protein", state="I")' % Y
            prod = '%s(gene="protein", state="A")' % Y
        if Y in compound_monomers:
            enzyme = '%s(bf=None)' % Y
            sub = '%s(gene="protein", state="I")' % X
            prod = '%s(gene="protein", state="A")' % X

        tmp1 = 'catalyze(%s, %s, %s, [%s, %s, %s] )\n' % (
            enzyme, sub, prod, rate1, rate2, rate3)

        rules += tmp1
        parameters += 'Parameter("%s", 1)\n' % rate1
        parameters += 'Parameter("%s", 1)\n' % rate2
        parameters += 'Parameter("%s", 1)\n' % rate3

    elif split[0] == 'binding/association':

        rate1 = "kf%s" % str(param_counter)
        rate2 = "kr%s" % str(param_counter)
        rule_name = '"%s_binds_%s_%s"' % (X, Y, str(param_counter))

        tmp1 = 'Rule(%s, %s(bf=None, gene="protein") + ' \
               '%s(bf=None, gene="protein") !> ' \
               '%s(bf=1, gene="protein") ** %s(bf=1, gene="protein"),' \
               ' %s, %s)\n' % (rule_name, X, Y, X, Y, rate1, rate2)

        tmp1 = tmp1.replace("**", "%")
        rules += tmp1
        parameters += 'Parameter("%s", 1)\n' % rate1
        parameters += 'Parameter("%s", 1)\n' % rate2

    else:
        print(split[0], 'not found')
    return rules, parameters


