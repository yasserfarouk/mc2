
# Repository: projects/mc2
echo -e "${BLUE}Restoring: projects/mc2${NC}"
if [ -d "projects/mc2/.git" ]; then
    echo -e "  ${YELLOW}Directory already exists, skipping...${NC}"
else
    # Create parent directory if needed
    mkdir -p "projects"
    
    # Clone the repository
    if git clone "git@github.com:yasserfarouk/mc2" "projects/mc2"; then
        echo -e "  ${GREEN}✓${NC} Successfully cloned"
        
        # Checkout the original branch if not already on it
        cd "projects/mc2"
        current=$(git rev-parse --abbrev-ref HEAD)
        if [ "$current" != "master" ]; then
            if git checkout "master" 2>/dev/null; then
                echo -e "  ${GREEN}✓${NC} Checked out branch: master"
            else
                echo -e "  ${YELLOW}⚠${NC} Could not checkout branch: master"
            fi
        fi
        cd - > /dev/null
    else
        echo -e "  ${RED}✗${NC} Failed to clone"
    fi
fi


# Repository: projects/mc2
echo -e "${BLUE}Restoring: projects/mc2${NC}"
if [ -d "projects/mc2/.git" ]; then
    echo -e "  ${YELLOW}Directory already exists, skipping...${NC}"
else
    # Create parent directory if needed
    mkdir -p "projects"
    
    # Clone the repository
    if git clone "git@github.com:yasserfarouk/mc2" "projects/mc2"; then
        echo -e "  ${GREEN}✓${NC} Successfully cloned"
        
        # Checkout the original branch if not already on it
        cd "projects/mc2"
        current=$(git rev-parse --abbrev-ref HEAD)
        if [ "$current" != "master" ]; then
            if git checkout "master" 2>/dev/null; then
                echo -e "  ${GREEN}✓${NC} Checked out branch: master"
            else
                echo -e "  ${YELLOW}⚠${NC} Could not checkout branch: master"
            fi
        fi
        cd - > /dev/null
    else
        echo -e "  ${RED}✗${NC} Failed to clone"
    fi
fi


# Repository: projects/mc2
echo -e "${BLUE}Restoring: projects/mc2${NC}"
if [ -d "projects/mc2/.git" ]; then
    echo -e "  ${YELLOW}Directory already exists, skipping...${NC}"
else
    # Create parent directory if needed
    mkdir -p "projects"
    
    # Clone the repository
    if git clone "git@github.com:yasserfarouk/mc2" "projects/mc2"; then
        echo -e "  ${GREEN}✓${NC} Successfully cloned"
        
        # Checkout the original branch if not already on it
        cd "projects/mc2"
        current=$(git rev-parse --abbrev-ref HEAD)
        if [ "$current" != "master" ]; then
            if git checkout "master" 2>/dev/null; then
                echo -e "  ${GREEN}✓${NC} Checked out branch: master"
            else
                echo -e "  ${YELLOW}⚠${NC} Could not checkout branch: master"
            fi
        fi
        cd - > /dev/null
    else
        echo -e "  ${RED}✗${NC} Failed to clone"
    fi
fi

