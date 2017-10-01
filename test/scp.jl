function wscope(init_x::T) where T
    need_x = true
    i = 0
    while i < 2
        if need_x
	    x = copy(init_x)
	    need_x = false
	end
        a = x[1]
        println(a)
	x[1] = x[1] + 2
	i += 1
    end

end

y = [1.,2,3,4]

wscope(y)
	