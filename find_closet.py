
def find_closet_number(data_list, target):
    """ 
    Finds the number in a list that is close3t  to  agiven value

    Args: 
    Data list:  a list of numbers. 
    target: the targuet number in to find the closet value to. 

    Retuns: 
    The number from data_list thet is closet to the target. 
    """

    if not data_list:
        return None
    index = min(range(len(data_list)), key=lambda i: abs(data_list[i] - target))
    return data_list[index], index
